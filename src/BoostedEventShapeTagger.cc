/*
Created:        28 November 2017
Last Updated:    4 December 2017

Justin Pilot
UC Davis

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University

-----

BoostedEventShapeTagger

Development of standalone BEST package for analyzers.
Based on the framework created by Justin:
  https://github.com/justinrpilot/BESTAnalysis

Requires MiniAOD inputs to access proper set of information

*/
#include "../interface/BoostedEventShapeTagger.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"
#include "PhysicsTools/CandUtils/interface/Thrust.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

// FASTJET
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/tools/Filter.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>


// constructor
BoostedEventShapeTagger::BoostedEventShapeTagger(const std::string& configFile, const std::string& dnnFile) :
  m_jetSoftDropMassMin(0),
  m_jetPtMin(0),
  m_radiusSmall(0),
  m_radiusLarge(0),
  m_reclusterJetPtMin(0),
  m_jetChargeKappa(0),
  m_maxJetSize(4){
    // Configuration
    std::vector<std::string> configurations;
    read_file(configFile,configurations);
    setConfigurations(configurations);


    // Basic options
    m_numSubjetsMin      = std::stoi(m_configurations.at("numSubjetsMin"));      // min. number of subjets
    m_numDaughtersMin    = std::stoi(m_configurations.at("numDaughtersMin"));    // min. daughters
    m_jetSoftDropMassMin = std::stof(m_configurations.at("jetSoftDropMassMin")); // min. soft drop mass
    m_jetPtMin           = std::stof(m_configurations.at("jetPtMin"));           // min. jet pT

    m_radiusSmall        = std::stof(m_configurations.at("radiusSmall"));
    m_radiusLarge        = std::stof(m_configurations.at("radiusLarge"));
    m_reclusterJetPtMin  = std::stof(m_configurations.at("reclusterJetPtMin"));

    m_jetChargeKappa     = std::stof(m_configurations.at("jetChargeKappa"));
    m_maxJetSize         = std::stoi(m_configurations.at("maxJetSize"));

    // DNN material lwtnn interface
//    std::string dnnFile = m_configurations.at("dnnFile");
     std::cout << "Using dnn file " << dnnFile << std::endl;
    std::ifstream input_cfg( dnnFile );                     // original: "data/BEST_mlp.json"
    lwt::JSONConfig cfg = lwt::parse_json( input_cfg );
    m_lwtnn = std::make_unique<lwt::LightweightNeuralNetwork>(cfg.inputs, cfg.layers, cfg.outputs);
} // end constructor


BoostedEventShapeTagger::~BoostedEventShapeTagger(){
}


const std::map<std::string,double>& BoostedEventShapeTagger::execute( const pat::Jet& jet ){
    /* Dan Guest's lightweight DNN framework */
    getJetValues(jet);                            // update m_BESTvars
    if (m_BESTvars.empty()){
      m_NNresults["dnn_qcd"] = 1;
      m_NNresults["dnn_top"] = 0;
      m_NNresults["dnn_higgs"] = 0;
      m_NNresults["dnn_z"] = 0;
      m_NNresults["dnn_w"] = 0;
      m_NNresults["dnn_b"] = 0;
    }else{
      m_NNresults = m_lwtnn->compute(m_BESTvars);
    }
    return m_NNresults;
}


void BoostedEventShapeTagger::getJetValues( const pat::Jet& jet ){
    /* Grab attributes from the jet and store them in map
       Jet requirements:
         pT > 500 GeV
         soft-drop mass > 40 GeV
         >=2 subjets
         >=2 daughters
           daughter->pt() >= 0.05
         m_reclusterJetPtMin = 30 GeV
         maxJets = 4
           sumP, sumPz
           pair-wise invariant mass
    */
    using namespace fastjet;
    typedef reco::Candidate::PolarLorentzVector fourv;

    // clear the map from the previous jet's values
    m_BESTvars.clear();

    // Access the subjets
    auto const& thisSubjets   = jet.subjets();
    unsigned int numDaughters = jet.numberOfDaughters();

    // Do some checks on the jets and print warnings to the user
    if (thisSubjets.size() < m_numSubjetsMin){
        edm::LogInfo("BEST") << " WARNING :: BEST : Number of subjets, " << thisSubjets.size() << ", is less than " << m_numSubjetsMin << std::endl;
        edm::LogInfo("BEST") << " WARNING :: BEST : -- BEST will NOT run.  Please check your jets! " << std::endl;
        return;
    }
    if (numDaughters < m_numDaughtersMin){
        edm::LogInfo("BEST") << " WARNING :: BEST : Number of daughters, " << numDaughters << ", is less than " << m_numDaughtersMin << std::endl;
        edm::LogInfo("BEST") << " WARNING :: BEST : -- BEST will NOT run.  Please check your jets! " << std::endl;
        return;
    }
    if ( jet.pt() < m_jetPtMin){
        edm::LogInfo("BEST") << " WARNING :: BEST : The jet pT " << jet.pt() << ", is less than " << m_jetPtMin << std::endl;
        edm::LogInfo("BEST") << " WARNING :: BEST : -- BEST will run, but the results can't be trusted! Please check your jets! " << std::endl;
    }
    if (jet.userFloat("ak8PFJetsPuppiSoftDropMass") < m_jetSoftDropMassMin){
        edm::LogInfo("BEST") << " WARNING :: BEST : The soft-drop mass " << jet.userFloat("ak8PFJetsPuppiSoftDropMass") << ", is less than " << m_jetSoftDropMassMin << std::endl;
        edm::LogInfo("BEST") << " WARNING :: BEST : -- BEST will run, but the results can't be trusted! Please check your jets! " << std::endl;
    }


    // b-tagging
    float btagValue1 = thisSubjets.at(0)->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    float btagValue2 = thisSubjets.at(1)->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

    // n-subjettiness
    float tau1 = jet.userFloat("NjettinessAK8Puppi:tau1");
    float tau2 = jet.userFloat("NjettinessAK8Puppi:tau2");
    float tau3 = jet.userFloat("NjettinessAK8Puppi:tau3");

    // BEST vars
    fourv thisJet = jet.polarP4();

    // Boost to rest frame
    TLorentzVector thisJetLV;
    TLorentzVector thisJetLV_W;
    TLorentzVector thisJetLV_Z;
    TLorentzVector thisJetLV_H;
    TLorentzVector thisJetLV_top;
    thisJetLV.SetPtEtaPhiM(thisJet.Pt(), thisJet.Eta(), thisJet.Phi(), thisJet.M() );
    thisJetLV_W.SetPtEtaPhiM(thisJet.Pt(),   thisJet.Eta(), thisJet.Phi(), m_Wmass ); // W mass [GeV]
    thisJetLV_Z.SetPtEtaPhiM(thisJet.Pt(),   thisJet.Eta(), thisJet.Phi(), m_Zmass ); // Z mass
    thisJetLV_H.SetPtEtaPhiM(thisJet.Pt(),   thisJet.Eta(), thisJet.Phi(), m_Hmass ); // Higgs mass
    thisJetLV_top.SetPtEtaPhiM(thisJet.Pt(), thisJet.Eta(), thisJet.Phi(), m_Tmass ); // Top mass

    std::vector<TLorentzVector> particles_jet;
    std::vector<TLorentzVector> particles_top;
    std::vector<TLorentzVector> particles_W;
    std::vector<TLorentzVector> particles_Z;
    std::vector<TLorentzVector> particles_H;

    std::vector<math::XYZVector> particles2_jet;
    std::vector<math::XYZVector> particles2_top;
    std::vector<math::XYZVector> particles2_W;
    std::vector<math::XYZVector> particles2_Z;
    std::vector<math::XYZVector> particles2_H;

    std::vector<reco::LeafCandidate> particles3_jet;
    std::vector<reco::LeafCandidate> particles3_top;
    std::vector<reco::LeafCandidate> particles3_W;
    std::vector<reco::LeafCandidate> particles3_Z;
    std::vector<reco::LeafCandidate> particles3_H;

    std::vector<fastjet::PseudoJet> topFJparticles;
    std::vector<fastjet::PseudoJet> ZFJparticles;
    std::vector<fastjet::PseudoJet> WFJparticles;
    std::vector<fastjet::PseudoJet> HFJparticles;
    std::vector<fastjet::PseudoJet> topFJparticles_noBoost;
    std::vector<fastjet::PseudoJet> jetFJparticles;
    std::vector<fastjet::PseudoJet> jetFJparticles_noBoost;
    std::vector<fastjet::PseudoJet> jetFJparticles_transformed;

    TVector3 transformedV;

    // Jet charge calculation (from daughters)
    float qxptsum(0.0);                             // jet charge
    float ptsum = pow(jet.pt(), m_jetChargeKappa);  // weighted jet pT

    // loop over daughters
    std::vector<const pat::PackedCandidate*> daughtersOfJet;   // store all daughters in one vector

    for (unsigned int i=0,size=jet.numberOfDaughters(); i<size; i++){
        daughtersOfJet.push_back( dynamic_cast<const pat::PackedCandidate*>(jet.daughter(i)) );
    }

    for(unsigned int i=0,size=daughtersOfJet.size(); i<size; i++){

        const auto* daughter = daughtersOfJet.at(i);

        if (daughter->pt() < 0.5) continue;

        float dau_px = daughter->px();
        float dau_py = daughter->py();
        float dau_pz = daughter->pz();
        float dau_e  = daughter->energy();

        TLorentzVector thisParticleLV_jet( dau_px,dau_py,dau_pz,dau_e );
        TLorentzVector thisParticleLV_top( dau_px,dau_py,dau_pz,dau_e );
        TLorentzVector thisParticleLV_W(   dau_px,dau_py,dau_pz,dau_e );
        TLorentzVector thisParticleLV_Z(   dau_px,dau_py,dau_pz,dau_e );
        TLorentzVector thisParticleLV_H(   dau_px,dau_py,dau_pz,dau_e );



        topFJparticles_noBoost.push_back( PseudoJet( thisParticleLV_top.X(), thisParticleLV_top.Y(), thisParticleLV_top.Z(), thisParticleLV_top.T() ) );
        jetFJparticles_noBoost.push_back( PseudoJet( thisParticleLV_jet.X(), thisParticleLV_jet.Y(), thisParticleLV_jet.Z(), thisParticleLV_jet.T() ) );

        TLorentzVector thisParticleLV_transformed( transformedV.X(), transformedV.Y(), transformedV.Z(), thisParticleLV_jet.E() );

        jetFJparticles_transformed.push_back( PseudoJet( thisParticleLV_transformed.X(), thisParticleLV_transformed.Y(), thisParticleLV_transformed.Z(), thisParticleLV_transformed.T() ) );


        if (daughter->pt() > 1.0)
            qxptsum += daughter->charge() * pow( daughter->pt(), m_jetChargeKappa);


        thisParticleLV_jet.Boost( -thisJetLV.BoostVector() );
        thisParticleLV_Z.Boost(   -thisJetLV_Z.BoostVector() );
        thisParticleLV_W.Boost(   -thisJetLV_W.BoostVector() );
        thisParticleLV_H.Boost(   -thisJetLV_H.BoostVector() );
        thisParticleLV_top.Boost( -thisJetLV_top.BoostVector() );

        pboost( thisJetLV_W.Vect(),   thisParticleLV_W.Vect(), thisParticleLV_W);
        pboost( thisJetLV_Z.Vect(),   thisParticleLV_Z.Vect(), thisParticleLV_Z);
        pboost( thisJetLV_top.Vect(), thisParticleLV_top.Vect(), thisParticleLV_top);
        pboost( thisJetLV_H.Vect(),   thisParticleLV_H.Vect(), thisParticleLV_H);

        particles_jet.push_back( thisParticleLV_jet );
        particles_top.push_back( thisParticleLV_top );
        particles_W.push_back(   thisParticleLV_W );
        particles_Z.push_back(   thisParticleLV_Z );
        particles_H.push_back(   thisParticleLV_H );

        topFJparticles.push_back( PseudoJet( thisParticleLV_top.X(), thisParticleLV_top.Y(), thisParticleLV_top.Z(), thisParticleLV_top.T() ) );
        WFJparticles.push_back(   PseudoJet( thisParticleLV_W.X(), thisParticleLV_W.Y(), thisParticleLV_W.Z(), thisParticleLV_W.T() ) );
        ZFJparticles.push_back(   PseudoJet( thisParticleLV_Z.X(), thisParticleLV_Z.Y(), thisParticleLV_Z.Z(), thisParticleLV_Z.T() ) );
        HFJparticles.push_back(   PseudoJet( thisParticleLV_H.X(), thisParticleLV_H.Y(), thisParticleLV_H.Z(), thisParticleLV_H.T() ) );
        jetFJparticles.push_back( PseudoJet( thisParticleLV_jet.X(), thisParticleLV_jet.Y(), thisParticleLV_jet.Z(), thisParticleLV_jet.T() ) );

        particles2_top.push_back( math::XYZVector( thisParticleLV_top.X(), thisParticleLV_top.Y(), thisParticleLV_top.Z() ));
        particles3_top.push_back( reco::LeafCandidate(+1, reco::Candidate::LorentzVector( thisParticleLV_top.X(), thisParticleLV_top.Y(), thisParticleLV_top.Z(), thisParticleLV_top.T()     ) ));
        particles2_W.push_back(   math::XYZVector( thisParticleLV_W.X(), thisParticleLV_W.Y(), thisParticleLV_W.Z() ));
        particles3_W.push_back(   reco::LeafCandidate(+1, reco::Candidate::LorentzVector( thisParticleLV_W.X(), thisParticleLV_W.Y(), thisParticleLV_W.Z(), thisParticleLV_W.T()     ) ));
        particles2_Z.push_back(   math::XYZVector( thisParticleLV_Z.X(), thisParticleLV_Z.Y(), thisParticleLV_Z.Z() ));
        particles3_Z.push_back(   reco::LeafCandidate(+1, reco::Candidate::LorentzVector( thisParticleLV_Z.X(), thisParticleLV_Z.Y(), thisParticleLV_Z.Z(), thisParticleLV_Z.T()     ) ));
        particles2_H.push_back(   math::XYZVector( thisParticleLV_H.X(), thisParticleLV_H.Y(), thisParticleLV_H.Z() ));
        particles3_H.push_back(   reco::LeafCandidate(+1, reco::Candidate::LorentzVector( thisParticleLV_H.X(), thisParticleLV_H.Y(), thisParticleLV_H.Z(), thisParticleLV_H.T()     ) ));
    } // end loop over daughters

    float jetq = qxptsum / ptsum;  // Jet Charge

    // Fox-Wolfram Moments
    double fwm[5]     = {0.0, 0.0 ,0.0 ,0.0, 0.0};
    double fwm_W[5]   = {0.0, 0.0 ,0.0 ,0.0, 0.0};
    double fwm_top[5] = {0.0, 0.0 ,0.0 ,0.0, 0.0};
    double fwm_Z[5]   = {0.0, 0.0 ,0.0 ,0.0, 0.0};
    double fwm_H[5]   = {0.0, 0.0 ,0.0 ,0.0, 0.0};

    FWMoments( particles_W,   fwm_W);
    FWMoments( particles_jet, fwm);
    FWMoments( particles_top, fwm_top);
    FWMoments( particles_Z,   fwm_Z);
    FWMoments( particles_H,   fwm_H);

    // Event Shapes
    EventShapeVariables eventShapes_top( particles2_top );
    EventShapeVariables eventShapes_W( particles2_W );
    EventShapeVariables eventShapes_Z( particles2_Z );
    EventShapeVariables eventShapes_H( particles2_H );

    // Thrust
    Thrust thrustCalculator_top( particles3_top.begin(), particles3_top.end() );
    Thrust thrustCalculator_W( particles3_W.begin(), particles3_W.end() );
    Thrust thrustCalculator_Z( particles3_Z.begin(), particles3_Z.end() );
    Thrust thrustCalculator_H( particles3_H.begin(), particles3_H.end() );

    // Recluster constituents
    JetDefinition jet_def(antikt_algorithm,  m_radiusSmall);
    JetDefinition jet_def2(antikt_algorithm, m_radiusLarge);

    ClusterSequence cs(    topFJparticles, jet_def);
    ClusterSequence cs_W(  WFJparticles,   jet_def);
    ClusterSequence cs_Z(  ZFJparticles,   jet_def);
    ClusterSequence cs_H(  HFJparticles,   jet_def);
    ClusterSequence cs_jet(jetFJparticles, jet_def);
    ClusterSequence cs_noBoost(topFJparticles_noBoost, jet_def2);

    ClusterSequence cs_transformed(jetFJparticles_transformed, jet_def);

    std::vector<PseudoJet> jetsFJ     = sorted_by_pt( cs.inclusive_jets(m_reclusterJetPtMin) );
    std::vector<PseudoJet> jetsFJ_W   = sorted_by_pt( cs_W.inclusive_jets(m_reclusterJetPtMin) );
    std::vector<PseudoJet> jetsFJ_Z   = sorted_by_pt( cs_Z.inclusive_jets(m_reclusterJetPtMin) );
    std::vector<PseudoJet> jetsFJ_H   = sorted_by_pt( cs_H.inclusive_jets(m_reclusterJetPtMin) );
    std::vector<PseudoJet> jetsFJ_jet = sorted_by_pt( cs_jet.inclusive_jets(m_reclusterJetPtMin) );
    std::vector<PseudoJet> jetsFJ_noBoost     = sorted_by_pt( cs_noBoost.inclusive_jets(m_reclusterJetPtMin) );
    std::vector<PseudoJet> jetsFJ_transformed = sorted_by_pt( cs_transformed.inclusive_jets(m_reclusterJetPtMin) );


    // pair-wise invariant masses
    TLorentzVector m1234LV_jet(0.,0.,0.,0.);
    TLorentzVector m1234LV_W(0.,0.,0.,0.);
    TLorentzVector m1234LV_Z(0.,0.,0.,0.);
    TLorentzVector m1234LV_H(0.,0.,0.,0.);
    TLorentzVector m1234LV_top(0.,0.,0.,0.);

    TLorentzVector m12LV_jet(0.,0.,0.,0.);
    TLorentzVector m12LV_W(0.,0.,0.,0.);
    TLorentzVector m12LV_Z(0.,0.,0.,0.);
    TLorentzVector m12LV_H(0.,0.,0.,0.);
    TLorentzVector m12LV_top(0.,0.,0.,0.);

    TLorentzVector m13LV_jet(0.,0.,0.,0.);
    TLorentzVector m13LV_W(0.,0.,0.,0.);
    TLorentzVector m13LV_Z(0.,0.,0.,0.);
    TLorentzVector m13LV_H(0.,0.,0.,0.);
    TLorentzVector m13LV_top(0.,0.,0.,0.);

    TLorentzVector m23LV_jet(0.,0.,0.,0.);
    TLorentzVector m23LV_W(0.,0.,0.,0.);
    TLorentzVector m23LV_Z(0.,0.,0.,0.);
    TLorentzVector m23LV_H(0.,0.,0.,0.);
    TLorentzVector m23LV_top(0.,0.,0.,0.);

    // sum of jet pz and p  Indices = top, W, Z, H, j
    float sumP[5]  = {0.0,0.0,0.0,0.0,0.0};
    float sumPz[5] = {0.0,0.0,0.0,0.0,0.0};

    // -- top
    for (size_t ii=0,size=std::min(m_maxJetSize,jetsFJ.size()); ii<size; ii++){
        sumPz[0] += jetsFJ[ii].pz();
        sumP[0]  += sqrt( jetsFJ[ii].modp2() );
        thisJetLV    = TLorentzVector(jetsFJ[ii].px(), jetsFJ[ii].py(), jetsFJ[ii].pz(), jetsFJ[ii].e());
        m1234LV_top += thisJetLV;
        switch (ii){
            case 0:
                m12LV_top += thisJetLV;
                m13LV_top += thisJetLV;
                break;
            case 1:
                m12LV_top += thisJetLV;
                m23LV_top += thisJetLV;
                break;
            case 2:
                m13LV_top += thisJetLV;
                m23LV_top += thisJetLV;
                break;
            case 3:
                break;
        }
    }

    // -- W jets
    for (size_t ii=0,size=std::min(m_maxJetSize,jetsFJ_W.size()); ii<size; ii++){
        sumPz[1] += jetsFJ_W[ii].pz();
        sumP[1]  += sqrt( jetsFJ_W[ii].modp2() );
        thisJetLV  = TLorentzVector(jetsFJ_W[ii].px(), jetsFJ_W[ii].py(), jetsFJ_W[ii].pz(), jetsFJ_W[ii].e());
        m1234LV_W += thisJetLV;
        switch (ii){
            case 0:
                m12LV_W += thisJetLV;
                m13LV_W += thisJetLV;
                break;
            case 1:
                m12LV_W += thisJetLV;
                m23LV_W += thisJetLV;
                break;
            case 2:
                m13LV_W += thisJetLV;
                m23LV_W += thisJetLV;
                break;
            case 3:
                break;
        }
    }

    // -- Z jets
    for (size_t ii=0,size=std::min(m_maxJetSize,jetsFJ_Z.size()); ii<size; ii++){
        sumPz[2] += jetsFJ_Z[ii].pz();
        sumP[2]  += sqrt( jetsFJ_Z[ii].modp2() );
        thisJetLV  = TLorentzVector(jetsFJ_Z[ii].px(), jetsFJ_Z[ii].py(), jetsFJ_Z[ii].pz(), jetsFJ_Z[ii].e());
        m1234LV_Z += thisJetLV;
        switch (ii){
            case 0:
                m12LV_Z += thisJetLV;
                m13LV_Z += thisJetLV;
                break;
            case 1:
                m12LV_Z += thisJetLV;
                m23LV_Z += thisJetLV;
                break;
            case 2:
                m13LV_Z += thisJetLV;
                m23LV_Z += thisJetLV;
                break;
            case 3:
                break;
        }
    }

    // -- Higgs jets
    for (size_t ii=0,size=std::min(m_maxJetSize,jetsFJ_H.size()); ii<size; ii++){
        sumPz[3] += jetsFJ_H[ii].pz();
        sumP[3]  += sqrt( jetsFJ_H[ii].modp2() );
        thisJetLV  = TLorentzVector(jetsFJ_H[ii].px(), jetsFJ_H[ii].py(), jetsFJ_H[ii].pz(), jetsFJ_H[ii].e());
        m1234LV_H += thisJetLV;
        switch (ii){
            case 0:
                m12LV_H += thisJetLV;
                m13LV_H += thisJetLV;
                break;
            case 1:
                m12LV_H += thisJetLV;
                m23LV_H += thisJetLV;
                break;
            case 2:
                m13LV_H += thisJetLV;
                m23LV_H += thisJetLV;
                break;
            case 3:
                break;
        }
    }

    // -- jets
    for (size_t ii=0,size=std::min(m_maxJetSize,jetsFJ_jet.size()); ii<size; ii++){
        sumPz[4] += jetsFJ_jet[ii].pz();
        sumP[4]  += sqrt( jetsFJ_jet[ii].modp2() );
        thisJetLV    = TLorentzVector(jetsFJ_jet[ii].px(), jetsFJ_jet[ii].py(), jetsFJ_jet[ii].pz(), jetsFJ_jet[ii].e());
        m1234LV_jet += thisJetLV;
        switch (ii){
            case 0:
                m12LV_jet += thisJetLV;
                m13LV_jet += thisJetLV;
                break;
            case 1:
                m12LV_jet += thisJetLV;
                m23LV_jet += thisJetLV;
                break;
            case 2:
                m13LV_jet += thisJetLV;
                m23LV_jet += thisJetLV;
                break;
            case 3:
                break;
        }
    }


    // Update the map with new values
    m_BESTvars["bDisc"]   = (btagValue1 > btagValue2) ? btagValue1 : btagValue2;
    m_BESTvars["bDisc1"]  = btagValue1;
    m_BESTvars["bDisc2"]  = btagValue2;
    m_BESTvars["et"]      = thisJet.Pt();
    m_BESTvars["eta"]     = thisJet.Rapidity();
    m_BESTvars["mass"]    = thisJet.M();
    m_BESTvars["SDmass"]  = jet.groomedMass();
    m_BESTvars["tau32"]   = (tau2 > 1e-8) ? tau3/tau2 : 999.;
    m_BESTvars["tau21"]   = (tau1 > 1e-8) ? tau2/tau1 : 999.;
    m_BESTvars["q"]       = jetq;

    m_BESTvars["m1234_jet"] = m1234LV_jet.M();
    m_BESTvars["m12_jet"]   = m12LV_jet.M();
    m_BESTvars["m23_jet"]   = m23LV_jet.M();
    m_BESTvars["m13_jet"]   = m13LV_jet.M();

    m_BESTvars["m1234_top"] = m1234LV_top.M();
    m_BESTvars["m12_top"]   = m12LV_top.M();
    m_BESTvars["m23_top"]   = m23LV_top.M();
    m_BESTvars["m13_top"]   = m13LV_top.M();

    m_BESTvars["m1234_W"] = m1234LV_W.M();
    m_BESTvars["m12_W"]   = m12LV_W.M();
    m_BESTvars["m23_W"]   = m23LV_W.M();
    m_BESTvars["m13_W"]   = m13LV_W.M();

    m_BESTvars["m1234_Z"] = m1234LV_Z.M();
    m_BESTvars["m12_Z"]   = m12LV_Z.M();
    m_BESTvars["m23_Z"]   = m23LV_Z.M();
    m_BESTvars["m13_Z"]   = m13LV_Z.M();

    m_BESTvars["m1234_H"] = m1234LV_H.M();
    m_BESTvars["m12_H"]   = m12LV_H.M();
    m_BESTvars["m23_H"]   = m23LV_H.M();
    m_BESTvars["m13_H"]   = m13LV_H.M();

    std::vector<std::string> jetNames = {"top","W","Z","H","jet"};

    for (unsigned int pp=0; pp<5; pp++){
        std::string jetName = jetNames[pp];

        m_BESTvars["sumPz_"+jetName]   = sumPz[pp];
        m_BESTvars["sumP_"+jetName]    = sumP[pp];
        m_BESTvars["pzOverp_"+jetName] =  ( sumPz[pp] / (sumP[pp] + 0.0001) ); // not used for 'jet'
    }

    m_BESTvars["Njets_top"]  = jetsFJ.size();
    m_BESTvars["Njets_W"]    = jetsFJ_W.size();
    m_BESTvars["Njets_Z"]    = jetsFJ_Z.size();
    m_BESTvars["Njets_H"]    = jetsFJ_H.size();
    m_BESTvars["Njets_jet"]  = jetsFJ_jet.size();
    m_BESTvars["Njets_orig"] = jetsFJ_noBoost.size();

    // -- top values
    m_BESTvars["h1_top"] = fwm_top[1];
    m_BESTvars["h2_top"] = fwm_top[2];
    m_BESTvars["h3_top"] = fwm_top[3];
    m_BESTvars["h4_top"] = fwm_top[4];
    m_BESTvars["isotropy_top"]   = eventShapes_top.isotropy();
    m_BESTvars["sphericity_top"] = eventShapes_top.sphericity();
    m_BESTvars["aplanarity_top"] = eventShapes_top.aplanarity();
    m_BESTvars["thrust_top"]     = thrustCalculator_top.thrust();

    // -- W values
    m_BESTvars["h1_W"] = fwm_W[1];
    m_BESTvars["h2_W"] = fwm_W[2];
    m_BESTvars["h3_W"] = fwm_W[3];
    m_BESTvars["h4_W"] = fwm_W[4];
    m_BESTvars["isotropy_W"]   = eventShapes_W.isotropy();
    m_BESTvars["sphericity_W"] = eventShapes_W.sphericity();
    m_BESTvars["aplanarity_W"] = eventShapes_W.aplanarity();
    m_BESTvars["thrust_W"]     = thrustCalculator_W.thrust();

    // -- Z values
    m_BESTvars["h1_Z"] = fwm_Z[1];
    m_BESTvars["h2_Z"] = fwm_Z[2];
    m_BESTvars["h3_Z"] = fwm_Z[3];
    m_BESTvars["h4_Z"] = fwm_Z[4];
    m_BESTvars["isotropy_Z"]   = eventShapes_Z.isotropy();
    m_BESTvars["sphericity_Z"] = eventShapes_Z.sphericity();
    m_BESTvars["aplanarity_Z"] = eventShapes_Z.aplanarity();
    m_BESTvars["thrust_Z"]     = thrustCalculator_Z.thrust();

    // -- H values
    m_BESTvars["h1_H"] = fwm_H[1];
    m_BESTvars["h2_H"] = fwm_H[2];
    m_BESTvars["h3_H"] = fwm_H[3];
    m_BESTvars["h4_H"] = fwm_H[4];
    m_BESTvars["isotropy_H"]   = eventShapes_H.isotropy();
    m_BESTvars["sphericity_H"] = eventShapes_H.sphericity();
    m_BESTvars["aplanarity_H"] = eventShapes_H.aplanarity();
    m_BESTvars["thrust_H"]     = thrustCalculator_H.thrust();

    return;
}


void BoostedEventShapeTagger::pboost( TVector3 pbeam, TVector3 plab, TLorentzVector &pboo ){
    /* Given jet constituent momentum plab, find momentum relative to
       beam direction pbeam
    */
    double pl = plab.Dot(pbeam);
    pl *= 1 / pbeam.Mag();
    // double pt = sqrt(plab.Mag()*plab.Mag()-pl*pl);

    // set x axis direction along pbeam x (0,0,1)
    TVector3 pbx;
    pbx.SetX(pbeam.Y());
    pbx.SetY(pbeam.X());
    pbx.SetZ(0.0);
    pbx *= (1/pbx.Mag());

    // set y axis direction along -pbx x pbeam
    TVector3 pby;
    pby  = -pbx.Cross(pbeam);
    pby *= (1/pby.Mag());

    pboo.SetX(plab.Dot(pbx));
    pboo.SetY(plab.Dot(pby));
    pboo.SetZ(pl);

    return;
}


void BoostedEventShapeTagger::FWMoments( std::vector<TLorentzVector> particles, double (&outputs)[5] ){
    /* Fox-Wolfram moments */
    int numParticles = particles.size();

    float s(0.0);
    for(int i = 0; i < numParticles; i++){
        s += particles[i].E();
    }

    float H0(0.0);
    float H4(0.0);
    float H3(0.0);
    float H2(0.0);
    float H1(0.0);

    for (int i=0; i<numParticles; i++){
        for (int j=i; j<numParticles; j++){
            float costh = ( particles[i].Px() * particles[j].Px() + particles[i].Py() * particles[j].Py() + particles[i].Pz() * particles[j].Pz() ) / ( particles[i].P() * particles[j].P() );
            float w1 = particles[i].P();
            float w2 = particles[j].P();

            float fw0 = LegP(costh, 0);
            float fw1 = LegP(costh, 1);
            float fw2 = LegP(costh, 2);
            float fw3 = LegP(costh, 3);
            float fw4 = LegP(costh, 4);

            H0 += w1 * w2 * fw0;
            H1 += w1 * w2 * fw1;
            H2 += w1 * w2 * fw2;
            H3 += w1 * w2 * fw3;
            H4 += w1 * w2 * fw4;
        }
    }

    H0 += 1e-3;              // prevent dividing by 0
    outputs[0] = (H0);
    outputs[1] = (H1 / H0);
    outputs[2] = (H2 / H0);
    outputs[3] = (H3 / H0);
    outputs[4] = (H4 / H0);

    return;
}


float BoostedEventShapeTagger::LegP(float x, int order){
    /* Calculation in FWMoments */
    float value(0.0);

    if (order == 0) value = 1;
    else if (order == 1) value = x;
    else if (order == 2) value = 0.5*(3*x*x - 1);
    else if (order == 3) value = 0.5*(5*x*x*x - 3*x);
    else if (order == 4) value = 0.125*(35*x*x*x*x - 30*x*x + 3);
    else value = 0;

    return value;
}


unsigned int BoostedEventShapeTagger::getParticleID(){
    /* Use simple algorithm to get the predicted particle ID
       - Particle ID = Particle Type with largest score
           (particleType == 0) QCD
           (particleType == 1) Top
           (particleType == 2) H
           (particleType == 3) Z
           (particleType == 4) W

        Here you can also add more sophisticated algorithms for determining the tagging,
        e.g., define working points to "tag" a jet.
    */
    std::vector<double> values{ m_NNresults["dnn_qcd"],   m_NNresults["dnn_top"],
                                m_NNresults["dnn_higgs"], m_NNresults["dnn_z"], m_NNresults["dnn_w"], m_NNresults["dnn_b"] };

    unsigned int particleID(0);
    double max_value(-1.0);
    for (unsigned int pid=0,size=values.size();pid<size;pid++){
        if (values.at(pid) > max_value){
            max_value  = values.at(pid);
            particleID = pid;
        }
    }

    return particleID;
}


void BoostedEventShapeTagger::setConfigurations(const std::vector<std::string>& configurations){
    /* Set the configuration options */
    m_configurations.clear();
    for (const auto& config : configurations){
        // split config items by space
        std::istringstream cfg(config);
        std::istream_iterator<std::string> start(cfg), stop;
        std::vector<std::string> tokens(start, stop);

        m_configurations.insert( std::pair<std::string,std::string>(tokens.at(0),tokens.at(1)) );
    }

    // Protection against default settings missing in custom configuration
    // -- map of defaultConfigs defined in header
    for (const auto& defaultConfig : m_defaultConfigs){
        if ( m_configurations.find(defaultConfig.first) == m_configurations.end() ){
            // item isn't in config file
            edm::LogInfo("BEST") << " WARNING :: BEST : Configuration " << defaultConfig.first << " not defined." << std::endl;
            edm::LogInfo("BEST") << " WARNING :: BEST : Setting value to default value: " << defaultConfig.second << std::endl;
            m_configurations[defaultConfig.first] = defaultConfig.second;
        }
    }

    return;
}


void BoostedEventShapeTagger::read_file( const std::string &file_name, std::vector<std::string> &values, const std::string &comment ) {
    /* Read in a generic file and put it into a vector of strings */
    std::ifstream tmp_name(file_name.c_str());

    // open the file and put the data into a vector
    std::string line("");
    if (tmp_name.is_open()){

        while (std::getline(tmp_name, line)) {
            std::string newstring(line);

            // allow for comments
            std::size_t lineComment = line.find(comment);
            if (lineComment != std::string::npos) newstring = line.substr(0, lineComment);

            // remove all white spaces at the end of the string
            std::size_t space_pos = newstring.rfind(" ");
            while ( space_pos != std::string::npos && space_pos == newstring.size()-1 ) {
                newstring = newstring.substr(0, newstring.rfind(" "));
                space_pos = newstring.rfind(" ");
            }

            // ignore empty lines
            if(newstring.length()==0) continue;

            values.push_back(newstring); // put values into vector
        }

        tmp_name.close();
    }

    return;
}


bool BoostedEventShapeTagger::str2bool( const std::string value ){
    /* Turn string into boolean */
    bool valueBoolean(false);

    if (value.compare("True")==0 || value.compare("true")==0 || value.compare("1")==0){
        valueBoolean = true;
    }
    else{
        valueBoolean = false;
    }

    return valueBoolean;
}


// THE END
