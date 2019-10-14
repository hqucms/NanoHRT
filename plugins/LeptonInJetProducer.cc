// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "fastjet/PseudoJet.hh"
#include <fastjet/JetDefinition.hh>
#include <TLorentzVector.h>
#include <TMath.h>

#include "PhysicsTools/NanoAOD/interface/MatchingUtils.h"

// return index of lepton collection for hardest lepton in jet
// either by looking into the PF constituents                                                                                                                                  
// or by matching with dR to jet and constituents 
template <class C1>
void leptonInJet(const pat::Jet &jet, const C1 & leptons, float &lepdR, float &lepPt, int &lepIndex, int &lepId, float jetdR=0.8){
  float tmpdR=lepdR,tmpPt=lepPt;
  int tmpIndex(-1),tmpId(-1);

  for(unsigned ilep(0); ilep < leptons->size(); ilep++){
    auto itLep = leptons->ptrAt(ilep);
    float dR = reco::deltaR(jet.eta(), jet.phi(), itLep->eta(), itLep->phi());
    if( dR < tmpdR && dR < jetdR && itLep->pt() > tmpPt) {
      tmpdR = dR;
      tmpIndex = ilep;
      tmpPt = itLep->pt();
      if(itLep->isMuon()) tmpId = 13;
      if(itLep->isElectron()) tmpId = 11;
      break;
    }
  } // loop over leptons

  bool matched = false;
  if(tmpIndex>-1){
    auto itLep = leptons->ptrAt(tmpIndex);
    if(matchByCommonSourceCandidatePtr(*itLep,jet)) { matched =true; lepId = tmpId; std::cout << " matched by source " << std::endl;}
    else{
      for (auto const pPart : jet.daughterPtrVector() ) {
	if(reco::deltaR(pPart->eta(), pPart->phi(), itLep->eta(), itLep->phi()) < 0.01){
	  if(fabs(pPart->pdgId()) == tmpId){
	    matched =true; lepId = tmpId;
	    break;
	  }
	  else if(fabs(pPart->pt() - itLep->pt()) < 0.01){
	    matched =true; lepId = fabs(pPart->pdgId()); // for now even save those that are PFeles and muons
	    break;
	  }
	}
      } // loop over jet daughters
    }
    if(matched){
      lepdR = tmpdR;
      lepPt = tmpPt;
      lepIndex = tmpIndex;
    }
  }// check if matched

}

template <typename T>
class LeptonInJetProducer : public edm::global::EDProducer<> {
public:
  explicit LeptonInJetProducer(const edm::ParameterSet &iConfig) :
    srcJet_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("src"))),
    srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
    srcMu_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("srcMu")))
  {
    produces<edm::ValueMap<float>>("lsf3");
    produces<edm::ValueMap<float>>("dRLep");
    produces<edm::ValueMap<int>>("muIdx3SJ");
    produces<edm::ValueMap<int>>("eleIdx3SJ");
    produces<edm::ValueMap<int>>("idLep");
  }
  ~LeptonInJetProducer() override {};
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override ;

  static bool orderPseudoJet(fastjet::PseudoJet j1, fastjet::PseudoJet j2);
  std::tuple<float,float> calculateLSF(std::vector<fastjet::PseudoJet> iCParticles, std::vector<fastjet::PseudoJet> &ljets,
				       float ilPt, float ilEta, float ilPhi, int ilId, double dr, int nsj) const;

  edm::EDGetTokenT<edm::View<pat::Jet>> srcJet_;
  edm::EDGetTokenT<edm::View<pat::Electron>> srcEle_;
  edm::EDGetTokenT<edm::View<pat::Muon>> srcMu_;

};

// ------------ method called to produce the data  ------------
template <typename T>
void
LeptonInJetProducer<T>::produce(edm::StreamID streamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const
{

    edm::Handle<edm::View<pat::Jet>> srcJet;
    iEvent.getByToken(srcJet_, srcJet);
    edm::Handle<edm::View<pat::Electron>> srcEle;
    iEvent.getByToken(srcEle_, srcEle);   
    edm::Handle<edm::View<pat::Muon>> srcMu;
    iEvent.getByToken(srcMu_, srcMu);

    std::vector<float> *vlsf3 = new std::vector<float>;
    std::vector<float> *vdRLep = new std::vector<float>;
    std::vector<int> *vmuIdx3SJ = new std::vector<int>;
    std::vector<int> *veleIdx3SJ = new std::vector<int>;
    std::vector<int> *vidLep = new std::vector<int>;

    // Find leptons in jets
    for (unsigned int ij = 0; ij<srcJet->size(); ij++){
      const pat::Jet &itJet = (*srcJet)[ij];
      if(itJet.pt() <= 10){
        vlsf3->push_back( -1);
        vdRLep->push_back( -1);
	veleIdx3SJ->push_back( -1);
        vmuIdx3SJ->push_back( -1);
	vidLep->push_back( -1);
        continue;
      }

      // take the jet constituents as the particles to compute the lsf
      std::vector<fastjet::PseudoJet> lClusterParticles;
      for (auto const d : itJet.daughterPtrVector() ) {
	fastjet::PseudoJet p( d->px(), d->py(), d->pz(), d->energy() );
        lClusterParticles.emplace_back(p);
      }

      // find the hard lepton inside the jet
      int ele_pfmatch_index(-1),mu_pfmatch_index(-1);
      float lepPt(-1),lepEta(-1),lepPhi(-1); 
      int lepId(-1); // save the id of the pf part matched for now
      float dRmin(0.8),dRtmp(999.),ptmin(10),dRLep(-1); //ptmin of lepton is 10 GeV

      leptonInJet(itJet, srcEle, dRtmp, ptmin, ele_pfmatch_index, lepId, dRmin);
      if(ele_pfmatch_index!=-1) {
        auto itLep = srcEle->ptrAt(ele_pfmatch_index);
	lepPt = itLep->pt(); lepEta = itLep->eta(); lepPhi = itLep->phi(); lepId = 11;
      }

      leptonInJet(itJet, srcMu, dRtmp, ptmin, mu_pfmatch_index, lepId, dRmin);
      if(mu_pfmatch_index!=-1) {
        auto itLep = srcMu->ptrAt(mu_pfmatch_index);
	lepPt =itLep->pt(); lepEta = itLep->eta(); lepPhi = itLep->phi(); lepId = 13;
      }

      std::vector<fastjet::PseudoJet> psub_3;
      std::sort(lClusterParticles.begin(),lClusterParticles.end(),orderPseudoJet);
      auto lsf_3 = calculateLSF(lClusterParticles, psub_3, lepPt, lepEta, lepPhi, lepId, 2.0, 3);
      vlsf3->push_back( std::get<0>(lsf_3));
      vdRLep->push_back( dRLep );
      veleIdx3SJ->push_back( ele_pfmatch_index );
      vmuIdx3SJ->push_back( mu_pfmatch_index );
      vidLep->push_back( lepId );
    }


    // Filling table
    std::unique_ptr<edm::ValueMap<float>> lsf3V(new edm::ValueMap<float>());
    edm::ValueMap<float>::Filler fillerlsf3(*lsf3V);
    fillerlsf3.insert(srcJet,vlsf3->begin(),vlsf3->end());
    fillerlsf3.fill();
    iEvent.put(std::move(lsf3V),"lsf3");

    std::unique_ptr<edm::ValueMap<float>> dRLepV(new edm::ValueMap<float>());
    edm::ValueMap<float>::Filler fillerdRLep(*dRLepV);
    fillerdRLep.insert(srcJet,vdRLep->begin(),vdRLep->end());
    fillerdRLep.fill();
    iEvent.put(std::move(dRLepV),"dRLep");

    std::unique_ptr<edm::ValueMap<int>> muIdx3SJV(new edm::ValueMap<int>());
    edm::ValueMap<int>::Filler fillermuIdx3SJ(*muIdx3SJV);
    fillermuIdx3SJ.insert(srcJet,vmuIdx3SJ->begin(),vmuIdx3SJ->end());
    fillermuIdx3SJ.fill();
    iEvent.put(std::move(muIdx3SJV),"muIdx3SJ");

    std::unique_ptr<edm::ValueMap<int>> eleIdx3SJV(new edm::ValueMap<int>());
    edm::ValueMap<int>::Filler fillereleIdx3SJ(*eleIdx3SJV);
    fillereleIdx3SJ.insert(srcJet,veleIdx3SJ->begin(),veleIdx3SJ->end());
    fillereleIdx3SJ.fill();
    iEvent.put(std::move(eleIdx3SJV),"eleIdx3SJ");

    std::unique_ptr<edm::ValueMap<int>> idLepV(new edm::ValueMap<int>());
    edm::ValueMap<int>::Filler filleridLep(*idLepV);
    filleridLep.insert(srcJet,vidLep->begin(),vidLep->end());
    filleridLep.fill();
    iEvent.put(std::move(idLepV),"idLep");
}

template <typename T>
bool LeptonInJetProducer<T>::orderPseudoJet(fastjet::PseudoJet j1, fastjet::PseudoJet j2) {
  return j1.perp2() > j2.perp2();
}

template <typename T>
std::tuple<float,float> LeptonInJetProducer<T>::calculateLSF(std::vector<fastjet::PseudoJet> iCParticles, std::vector<fastjet::PseudoJet> &lsubjets, 
							     float ilPt, float ilEta, float ilPhi, int ilId, double dr, int nsj) const {
  float lsf(-1),lmd(-1);
  if(ilPt>0 && ilId >-1) {    
    TLorentzVector ilep; 
    if(ilId == 11) ilep.SetPtEtaPhiM(ilPt, ilEta, ilPhi, 0.000511);
    else if(ilId == 13) ilep.SetPtEtaPhiM(ilPt, ilEta, ilPhi, 0.105658);
    else ilep.SetPtEtaPhiM(ilPt, ilEta, ilPhi, 0); // allow for non-PF leptons for now..
    fastjet::JetDefinition lCJet_def(fastjet::kt_algorithm, dr);
    fastjet::ClusterSequence lCClust_seq(iCParticles, lCJet_def);
    if (dr > 0.5) {
      lsubjets = sorted_by_pt(lCClust_seq.exclusive_jets_up_to(nsj));
    }
    else {
      lsubjets = sorted_by_pt(lCClust_seq.inclusive_jets());
    }
    int lId(-1);
    double dRmin = 999.;
    for (unsigned int i0=0; i0<lsubjets.size(); i0++) {
      double dR = reco::deltaR(lsubjets[i0].eta(), lsubjets[i0].phi(), ilep.Eta(), ilep.Phi());
      if ( dR < dRmin ) {
	dRmin = dR;
	lId = i0;
	}
    }
    if(lId != -1) {
      TLorentzVector pVec; pVec.SetPtEtaPhiM(lsubjets[lId].pt(), lsubjets[lId].eta(), lsubjets[lId].phi(), lsubjets[lId].m());
      lsf = ilep.Pt()/pVec.Pt();
      lmd = (ilep-pVec).M()/pVec.M();
    }
  }
  return std::tuple<float,float>(lsf,lmd);
}

template <typename T>
void LeptonInJetProducer<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src")->setComment("jet input collection");
  desc.add<edm::InputTag>("srcEle")->setComment("electron input collection");
  desc.add<edm::InputTag>("srcMu")->setComment("muon input collection");
  std::string modname;
  modname+="LepIn";
  if (typeid(T) == typeid(pat::Jet)) modname+="Jet";
  modname+="Producer";
  descriptions.add(modname,desc);
}


typedef LeptonInJetProducer<pat::Jet> LepInJetProducer;

//define this as a plug-in
DEFINE_FWK_MODULE(LepInJetProducer);
