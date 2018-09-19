#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

#include "NNKit/FatJetNN/interface/FatJetNN.h"

#include <memory>
#include <iostream>

using namespace deepntuples;

class DeepBoostedJetProducer : public edm::stream::EDProducer<> {

  public:
    explicit DeepBoostedJetProducer(const edm::ParameterSet&);
    ~DeepBoostedJetProducer() override;

    static void fillDescriptions(edm::ConfigurationDescriptions&);

  private:
    void beginStream(edm::StreamID) override {}
    void produce(edm::Event&, const edm::EventSetup&) override;
    void endStream() override {}

    const edm::EDGetTokenT<edm::View<pat::Jet>> src_;
    const bool has_puppi_weighted_daughters_;
    const double jet_radius_;
    std::string nominal_nn_path_;
    std::string decorr_nn_path_;

    std::unique_ptr<FatJetNN> fatjetNN_;
    std::unique_ptr<FatJetNN> decorrNN_;

};

DeepBoostedJetProducer::DeepBoostedJetProducer(const edm::ParameterSet& iConfig)
: src_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("src")))
, has_puppi_weighted_daughters_(iConfig.getParameter<bool>("hasPuppiWeightedDaughters"))
, jet_radius_(iConfig.getUntrackedParameter<double>("jet_radius", 0.8))
, nominal_nn_path_(iConfig.getUntrackedParameter<std::string>("nominal_nn_path", "NNKit/data/ak8/full"))
, decorr_nn_path_(iConfig.getUntrackedParameter<std::string>("decorrelated_nn_path", "NNKit/data/ak8/decorrelated"))
{

  // initialize the FatJetNN class in the constructor
  auto cc = consumesCollector();
  if (!nominal_nn_path_.empty()){
    fatjetNN_ = std::make_unique<FatJetNN>(iConfig, cc, jet_radius_);
    // load json for input variable transformation
    fatjetNN_->load_json(edm::FileInPath(nominal_nn_path_+"/preprocessing.json").fullPath());
    // load DNN model and parameter files
    fatjetNN_->load_model(edm::FileInPath(nominal_nn_path_+"/resnet-symbol.json").fullPath(),
        edm::FileInPath(nominal_nn_path_+"/resnet.params").fullPath());
  }

  if (!decorr_nn_path_.empty()){
    decorrNN_ = std::make_unique<FatJetNN>(iConfig, cc, jet_radius_);
    // load json for input variable transformation
    decorrNN_->load_json(edm::FileInPath(decorr_nn_path_+"/preprocessing.json").fullPath());
    // load DNN model and parameter files
    decorrNN_->load_model(edm::FileInPath(decorr_nn_path_+"/resnet-symbol.json").fullPath(),
        edm::FileInPath(decorr_nn_path_+"/resnet.params").fullPath());
  }

  produces<pat::JetCollection>();

}

DeepBoostedJetProducer::~DeepBoostedJetProducer(){
}

void DeepBoostedJetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void DeepBoostedJetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<edm::View<pat::Jet>> jets;
  iEvent.getByToken(src_, jets);

  if (fatjetNN_){
    fatjetNN_->readEvent(iEvent, iSetup);
  }
  if (decorrNN_){
    decorrNN_->readEvent(iEvent, iSetup);
  }

  auto outputs = std::make_unique<pat::JetCollection>();

  for (const auto &jet : *jets){
    JetHelper jet_helper(&jet, has_puppi_weighted_daughters_);

    pat::Jet newJet(jet);
    if (fatjetNN_){
      const auto& nn = fatjetNN_->predict(jet_helper);
      std::string prefix = "DeepBoostedJet:nn_";
      for (unsigned i=0; i<nn.size(); ++i){
        newJet.addUserFloat(prefix+FatJetNNHelper::categoryName(static_cast<FatJetNNHelper::FatJetNNCategory>(i)), nn.at(i));
      }
    }

    if (decorrNN_){
      const auto& md = decorrNN_->predict(jet_helper);
      std::string prefix = "DeepBoostedJet:decorr_nn_";
      for (unsigned i=0; i<md.size(); ++i){
        newJet.addUserFloat(prefix+FatJetNNHelper::categoryName(static_cast<FatJetNNHelper::FatJetNNCategory>(i)), md.at(i));
      }
    }

    outputs->push_back(newJet);

  }

  // put into the event
  iEvent.put(std::move(outputs));

}

//define this as a plug-in
DEFINE_FWK_MODULE(DeepBoostedJetProducer);
