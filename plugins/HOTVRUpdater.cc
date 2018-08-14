#include <memory>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Candidate/interface/CandAssociation.h"

#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"

//
// class declaration
//

class HOTVRUpdater : public edm::stream::EDProducer<> {
  public:
    explicit HOTVRUpdater(const edm::ParameterSet&);
    ~HOTVRUpdater();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginStream(edm::StreamID) override;
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endStream() override;

    // ----------member data ---------------------------
    const edm::EDGetTokenT<edm::View<pat::Jet>> src_;
    const edm::EDGetTokenT<edm::View<pat::Jet>> subjetsSrc_; // subjets w/ JEC applied
};

//
// constructors and destructor
//
HOTVRUpdater::HOTVRUpdater(const edm::ParameterSet& iConfig)
: src_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("src")))
, subjetsSrc_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("subjets")))
{
  // We make both the fat jets and subjets, and we must store them as separate collections
  produces<pat::JetCollection>();
}


HOTVRUpdater::~HOTVRUpdater()
{

  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
HOTVRUpdater::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<edm::View<pat::Jet>> jets;
  iEvent.getByToken(src_, jets);

  edm::Handle<edm::View<pat::Jet>> correctedSubjets;
  iEvent.getByToken(subjetsSrc_, correctedSubjets);

  auto jetCollection = std::make_unique<pat::JetCollection>();

  for (unsigned ijet=0; ijet<jets->size(); ++ijet){
    const auto &jet = jets->at(ijet);
    pat::Jet newJet;
    newJet.setP4(jet.p4());
    newJet.setJetArea(jet.jetArea());
    for (const auto &name : jet.userFloatNames()){
      newJet.addUserFloat(name, jet.userFloat(name));
    }

    std::vector<edm::Ptr<pat::Jet>> newSubjets;
    reco::Candidate::LorentzVector newP4Sum;
    for (const auto &subjet : jet.subjets()) {
      for (unsigned k=0; k<correctedSubjets->size(); ++k){
        // match to updatedSubjets by deltaR
        if (reco::deltaR2(*subjet, correctedSubjets->at(k)) < 1e-6){
          newP4Sum += correctedSubjets->at(k).p4();
          newSubjets.push_back(correctedSubjets->ptrAt(k));
          break;
        }
      }
    }

    if (newSubjets.size()!=jet.subjets().size()){
      throw cms::Exception("RuntimeError") << "Subjet number changed from " << jet.subjets().size() << " to " << newSubjets.size();
    }
    // sort subjets by pt
    std::sort(newSubjets.begin(), newSubjets.end(), [](const edm::Ptr<pat::Jet> &a, const edm::Ptr<pat::Jet> &b){ return a->pt() > b->pt(); });
    newJet.addSubjets(newSubjets);

    // udpate HOTVR jet with the sum of subjets p4
    newJet.setP4(newP4Sum);

    if (newSubjets.size() >= 3){
      newJet.addUserFloat("fpt",  newSubjets[0]->pt() / newJet.pt());
      newJet.addUserFloat("mmin", std::min({
        (newSubjets[0]->p4()+newSubjets[1]->p4()).mass(),
        (newSubjets[0]->p4()+newSubjets[2]->p4()).mass(),
        (newSubjets[1]->p4()+newSubjets[2]->p4()).mass(),
      }));
    } else {
      newJet.addUserFloat("fpt",  0);
      newJet.addUserFloat("mmin", 0);
    }

    jetCollection->push_back(newJet);
  }

  iEvent.put(std::move(jetCollection));
}


// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
HOTVRUpdater::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
HOTVRUpdater::endStream() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HOTVRUpdater::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HOTVRUpdater);
