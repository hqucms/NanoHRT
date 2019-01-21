#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

class MuSubProducer : public edm::stream::EDProducer<> {

  public:
    explicit MuSubProducer(const edm::ParameterSet&);
    ~MuSubProducer() override;
    static void fillDescriptions(edm::ConfigurationDescriptions&);
    static void globalEndJob();

  private:
    void beginStream(edm::StreamID) override {}
    void produce(edm::Event&, const edm::EventSetup&) override;
    void endStream() override {}

    const edm::EDGetTokenT<edm::View<pat::Muon>> src_;
    const edm::EDGetTokenT<edm::View<pat::PackedCandidate>> pfcs_;
    const edm::EDGetTokenT<reco::VertexCollection> vtxs_;
    const double ptmin_;

};

MuSubProducer::MuSubProducer(const edm::ParameterSet& iConfig)
: 
src_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("src"))),
pfcs_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("pfcs"))),
vtxs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vtxs"))),
ptmin_(iConfig.getParameter<double>("ptmin"))
{
  produces<std::vector<pat::PackedCandidate>>("");
}
MuSubProducer::~MuSubProducer(){
}

void MuSubProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void MuSubProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<edm::View<pat::Muon>>muons;
  edm::Handle<edm::View<pat::PackedCandidate>>PackedCandidates;
  edm::Handle<reco::VertexCollection> vtxs;

  std::vector<pat::PackedCandidate> allpf;

  iEvent.getByToken(vtxs_, vtxs);
  iEvent.getByToken(src_, muons);
  iEvent.getByToken(pfcs_, PackedCandidates);
  for(const auto &cand : *PackedCandidates)   
	{
	bool keepcand=true;
	for (const auto &curmuon : *muons)
		{
		if(vtxs->size()>0)
			{
			if((curmuon.pt()>ptmin_) and (curmuon.isTightMuon(vtxs->at(0))) and (curmuon.numberOfSourceCandidatePtrs()>0))
				{
				const auto mucand = curmuon.sourceCandidatePtr(0);
				if ((abs(mucand->px()-cand.px())<0.0000001) and (abs(mucand->py()-cand.py())<0.0000001) and (abs(mucand->pz()-cand.pz())<0.0000001)) 
					{
					keepcand=false;
					break;
					}		
				}
			}
  		}
	if(keepcand) allpf.push_back(cand);
	}

  //Add to jet userfloats
  auto outputs = std::make_unique<std::vector<pat::PackedCandidate>>(allpf);
  iEvent.put(std::move(outputs),"");
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuSubProducer);
