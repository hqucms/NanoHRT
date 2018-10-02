#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

#include "../interface/BoostedEventShapeTagger.h"

#include <memory>
#include <iostream>

class BESTProducer : public edm::stream::EDProducer<> {

  public:
    explicit BESTProducer(const edm::ParameterSet&);
    ~BESTProducer() override;

    static void fillDescriptions(edm::ConfigurationDescriptions&);

  private:
    void beginStream(edm::StreamID) override {}
    void produce(edm::Event&, const edm::EventSetup&) override;
    void endStream() override {}

    const edm::EDGetTokenT<edm::View<pat::Jet>> src_;
    edm::FileInPath cfg_path_;
    edm::FileInPath dnn_path_;

    std::unique_ptr<BoostedEventShapeTagger> tagger_;

};

BESTProducer::BESTProducer(const edm::ParameterSet& iConfig)
: src_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("src")))
, cfg_path_(iConfig.getParameter<edm::FileInPath>("config_path"))
, dnn_path_(iConfig.getParameter<edm::FileInPath>("dnn_path"))
{

  tagger_ = std::make_unique<BoostedEventShapeTagger>(cfg_path_.fullPath(), dnn_path_.fullPath());

  produces<pat::JetCollection>();
}

BESTProducer::~BESTProducer(){
}

void BESTProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void BESTProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<edm::View<pat::Jet>> jets;
  iEvent.getByToken(src_, jets);

  auto outputs = std::make_unique<pat::JetCollection>();

  for (const auto &jet : *jets){
    pat::Jet newJet(jet);
    const auto& nn_results = tagger_->execute(jet);
    for (const auto &p : nn_results){
      newJet.addUserFloat("BEST:"+p.first, p.second);
    }
    outputs->push_back(newJet);
  }

  // put into the event
  iEvent.put(std::move(outputs));

}

//define this as a plug-in
DEFINE_FWK_MODULE(BESTProducer);
