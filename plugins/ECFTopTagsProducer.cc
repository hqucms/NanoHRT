// based on https://github.com/sidnarayanan/ECFTest/blob/master/Producer/plugins/TopTagProducer.cc

#include <memory>
#include <vector>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/MeasureDefinition.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include "../interface/EnergyCorrelations.h"
#include "RecoJets/JetAlgorithms/interface/HEPTopTaggerWrapperV2.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

using namespace fastjet;

//
// class declaration
//

class ECFTopTagsProducer : public edm::stream::EDProducer<> {
public:
  explicit ECFTopTagsProducer(const edm::ParameterSet&);
  ~ECFTopTagsProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginStream(edm::StreamID) override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endStream() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<pat::Jet>> src_token_;
  edm::FileInPath bdt_path_;
  std::string bdt_name_;

  std::unique_ptr<fastjet::JetDefinition> jetDef{nullptr};
  std::unique_ptr<fastjet::AreaDefinition> areaDef{nullptr};
  std::unique_ptr<fastjet::GhostedAreaSpec> activeArea{nullptr};
  std::unique_ptr<fastjet::contrib::SoftDrop> softDrop{nullptr};

  std::unique_ptr<fastjet::contrib::Njettiness> tauN{nullptr};
  std::unique_ptr<ECFCalculator> ecfcalc{nullptr};
  std::unique_ptr<fastjet::HEPTopTaggerV2> htt{nullptr};

  std::unique_ptr<TMVA::Reader> bdt;

  // variables
  std::vector<float> ecfs;
  float tau32sd = 0;
  float htt_frec = 0;
  float htt_mass = 0;
  float spectator= 1;

  float get_ecf(int o, int n, int ib) {
    return std::get<ECFCalculator::ecfP>(ecfcalc->access(ECFCalculator::pos_type(o-1, n-1, ib)));
  }

  void reset() {
    tau32sd  = 0;
    htt_frec = 0;
    htt_mass = 0;
    for (auto& e : ecfs)
      e = 0;
  }

};

//
// constructors and destructor
//
ECFTopTagsProducer::ECFTopTagsProducer(const edm::ParameterSet& iConfig)
: src_token_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("src")))
, bdt_path_(iConfig.getUntrackedParameter<edm::FileInPath>("bdt_path", edm::FileInPath("PhysicsTools/NanoHRT/data/ECFTopTag/top_ecfbdt_v8_BDT.weights.xml")))
, bdt_name_(iConfig.getUntrackedParameter<std::string>("bdt_name", "v0"))
, ecfs(11, 0)
{
  // soft drop
  const int activeAreaRepeats = 1;
  const double ghostArea = 0.01;
  const double ghostEtaMax = 7.0;
  const double radius = 1.5;
  const double sdZcut = 0.15;
  const double sdBeta = 1.;

  jetDef.reset(new JetDefinition(cambridge_algorithm,
      radius));
  activeArea.reset(new GhostedAreaSpec(ghostEtaMax,activeAreaRepeats,ghostArea));
  areaDef.reset(new AreaDefinition(active_area_explicit_ghosts,*activeArea));
  softDrop.reset(new contrib::SoftDrop(sdBeta,sdZcut,radius));

  // substructure
  tauN = std::make_unique<contrib::Njettiness>(contrib::OnePass_KT_Axes(),
      contrib::NormalizedMeasure(1.,radius));


  ecfcalc.reset(new ECFCalculator());

  bool optimalR=true; bool doHTTQ=false;
  double minSJPt=0.; double minCandPt=0.;
  double sjmass=30.; double mucut=0.8;
  double filtR=0.3; int filtN=5;
  int mode=4; double minCandMass=0.;
  double maxCandMass=9999999.; double massRatioWidth=9999999.;
  double minM23Cut=0.; double minM13Cut=0.;
  double maxM13Cut=9999999.;  bool rejectMinR=false;
  htt.reset(new fastjet::HEPTopTaggerV2(optimalR,doHTTQ,
      minSJPt,minCandPt,
      sjmass,mucut,
      filtR,filtN,
      mode,minCandMass,
      maxCandMass,massRatioWidth,
      minM23Cut,minM13Cut,
      maxM13Cut,rejectMinR));

  // BDTs
  bdt = std::make_unique<TMVA::Reader>();
  bdt->AddVariable("ecfN_1_2_20/pow(ecfN_1_2_10,2.00)",&ecfs[0]);
  bdt->AddVariable("ecfN_1_3_40/ecfN_2_3_20"          ,&ecfs[1]);
  bdt->AddVariable("ecfN_3_3_10/pow(ecfN_1_3_40,0.75)",&ecfs[2]);
  bdt->AddVariable("ecfN_3_3_10/pow(ecfN_2_3_20,0.75)",&ecfs[3]);
  bdt->AddVariable("ecfN_3_3_20/pow(ecfN_3_3_40,0.50)",&ecfs[4]);
  bdt->AddVariable("ecfN_1_4_20/pow(ecfN_1_3_10,2.00)",&ecfs[5]);
  bdt->AddVariable("ecfN_1_4_40/pow(ecfN_1_3_20,2.00)",&ecfs[6]);
  bdt->AddVariable("ecfN_2_4_05/pow(ecfN_1_3_05,2.00)",&ecfs[7]);
  bdt->AddVariable("ecfN_2_4_10/pow(ecfN_1_3_10,2.00)",&ecfs[8]);
  bdt->AddVariable("ecfN_2_4_10/pow(ecfN_2_3_05,2.00)",&ecfs[9]);
  bdt->AddVariable("ecfN_2_4_20/pow(ecfN_1_3_20,2.00)",&ecfs[10]);
  bdt->AddVariable("tau32SD",&tau32sd);
  bdt->AddVariable("htt_frec",&htt_frec);
  bdt->AddSpectator("eventNumber",&spectator);
  bdt->AddSpectator("runNumber",&spectator);
  bdt->AddSpectator("pt",&spectator);
  bdt->AddSpectator("mSD",&spectator);
  bdt->BookMVA(bdt_name_, bdt_path_.fullPath());

  produces<pat::JetCollection>();

}


ECFTopTagsProducer::~ECFTopTagsProducer(){
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void ECFTopTagsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // read input collection
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByToken(src_token_, jets);

  auto outputs = std::make_unique<pat::JetCollection>();

  for (const auto& jet : *jets) {
    pat::Jet newJet(jet);

    reset();
    // convert the jet into pseudojet
    std::vector<PseudoJet> vjet;
    for (const auto &cand : jet.daughterPtrVector()) {
      // candidate's momentum is already scaled by Puppi
      if (cand->pt() < 0.01) continue;
      vjet.emplace_back(cand->px(), cand->py(), cand->pz(), cand->energy());
    }

    ClusterSequenceArea seq(vjet, *jetDef, *areaDef);
    auto alljets(sorted_by_pt(seq.inclusive_jets(0.1)));
    if (alljets.size() == 0) {
      newJet.addUserFloat("ecf_0", -1);
      newJet.addUserFloat("httMass", -1);
      newJet.addUserFloat("httFRec", -1);
      newJet.addUserFloat("tau32sd", -1);
      newJet.addUserFloat("ecfTopTagBDT", -1);
      outputs->push_back(newJet);
      continue;
    }

    // soft drop the jet
    auto& pj(alljets[0]);
    auto sd((*softDrop)(pj));
    auto sdConstituents(sorted_by_pt(sd.constituents()));

    // compute tauNs
    tau32sd = tauN->getTau(3, sdConstituents) / tauN->getTau(2, sdConstituents);

    // compute ECFs
    unsigned nFilter = std::min(100, (int)sdConstituents.size());
    std::vector<PseudoJet> sdConstsFiltered(sdConstituents.begin(), sdConstituents.begin() + nFilter);
    ecfcalc->calculate(sdConstsFiltered);
    ecfs[0]  = get_ecf(1,2,2) / pow(get_ecf(1,2,1),2.00);
    ecfs[1]  = get_ecf(1,3,3) / get_ecf(2,3,2);
    ecfs[2]  = get_ecf(3,3,1) / pow(get_ecf(1,3,3),0.75);
    ecfs[3]  = get_ecf(3,3,1) / pow(get_ecf(2,3,2),0.75);
    ecfs[4]  = get_ecf(3,3,2) / pow(get_ecf(3,3,3),0.50);
    ecfs[5]  = get_ecf(1,4,2) / pow(get_ecf(1,3,1),2.00);
    ecfs[6]  = get_ecf(1,4,3) / pow(get_ecf(1,3,2),2.00);
    ecfs[7]  = get_ecf(2,4,0) / pow(get_ecf(1,3,0),2.00);
    ecfs[8]  = get_ecf(2,4,1) / pow(get_ecf(1,3,1),2.00);
    ecfs[9]  = get_ecf(2,4,1) / pow(get_ecf(2,3,0),2.00);
    ecfs[10] = get_ecf(2,4,2) / pow(get_ecf(1,3,2),2.00);

    // compute HEPTopTag
    PseudoJet httJet(htt->result(pj));
    if (httJet != 0) {
      auto* s(static_cast<fastjet::HEPTopTaggerV2Structure*>(httJet.structure_non_const_ptr()));
      htt_mass = s->top_mass();
      htt_frec = s->fRec();
    } else {
      htt_mass = 0;
      htt_frec = 0;
    }

    newJet.addUserFloat("ecf_0", ecfs[0]);
    newJet.addUserFloat("httMass", htt_mass);
    newJet.addUserFloat("httFRec", htt_frec);
    newJet.addUserFloat("tau32sd", tau32sd);
    bool sane = get_ecf(2, 4, 2) > 0; // sanity check - if anything hits 0 or below, this will
    newJet.addUserFloat("ecfTopTagBDT", sane ? bdt->EvaluateMVA(bdt_name_): -1.2);
    outputs->push_back(newJet);
  }

  assert(outputs->size() == jets->size());
  // put into the event
  iEvent.put(std::move(outputs));
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void ECFTopTagsProducer::beginStream(edm::StreamID) {
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void ECFTopTagsProducer::endStream() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ECFTopTagsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ECFTopTagsProducer);
