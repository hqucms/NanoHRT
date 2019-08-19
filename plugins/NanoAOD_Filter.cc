#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "TLorentzVector.h"

#include <iostream>
#include <fstream>

#include <memory>
#include <iostream>
#include <Python.h>



class  NanoAOD_Filter : public edm::stream::EDFilter<> {
public:
  explicit NanoAOD_Filter( const edm::ParameterSet & );   
private:
  const edm::EDGetTokenT<edm::View<pat::Jet>> srcAK8_;
  const edm::EDGetTokenT<edm::View<pat::Jet>> srcAK4_;
  const edm::EDGetTokenT<edm::View<pat::Muon>> srcmu_;
  const edm::EDGetTokenT<edm::View<pat::Electron>> srcele_;
  bool filter( edm::Event &, const edm::EventSetup & );
  void beginRun(edm::Run const&, edm::EventSetup const&);
  void beginJob() ;
  void endJob() ;
 };


NanoAOD_Filter::NanoAOD_Filter(const edm::ParameterSet& iConfig):
//, srcAK8_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("srcAK8")))
 srcAK4_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("srcAK4")))
//, srcmu_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("srcmu")))
//, srcele_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcele")))
 {   

 }


bool NanoAOD_Filter::filter( edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //edm::Handle<edm::View<pat::Jet>> jetsAK8;
  //iEvent.getByToken(srcAK8_, jetsAK8);

  edm::Handle<edm::View<pat::Jet>> jetsAK4;
  iEvent.getByToken(srcAK4_, jetsAK4);

  //This section selects HT>800
  float totht=0.0;
  for (const auto &AK4pfjet : *jetsAK4)
	{
	if (AK4pfjet.pt()>30.0) totht+=AK4pfjet.pt();
	//std::cout<<"AK4,ht "<<AK4pfjet.pt()<<" "<<totht<<std::endl;
	}
  if(totht<800.0) return 0;

  //This section selects at least two AK8 jets w/pt > 200
  //std::cout<<"NAK "<<jetsAK8->size()<<std::endl;
  //if (jetsAK8->size()>1)
	//{	
	//if ((jetsAK8->at(0).pt()<200.0) or (jetsAK8->at(1).pt()<200.0))return 0;
	//std::cout<<"AK8_1,AK8_2 "<<jetsAK8->at(0).pt()<<" "<<jetsAK8->at(1).pt()<<std::endl;
	//}
  //else return 0;

  //std::cout<<"PASS"<<std::endl;
  
  //edm::Handle<edm::View<pat::Muon>> mus;
  //iEvent.getByToken(srcmu_, mus);


  //edm::Handle<edm::View<pat::Electron>> eles;
  //iEvent.getByToken(srcele_, eles);

  //for (const auto &mu : *mus)
	//{
	//std::cout<<"New MU "<<mu.pt() << " " << mu.isLooseMuon() <<std::endl;
	//if(mu.pt()>60 and mu.isLooseMuon())
	//{
	//std::cout<<"FAIL mu pt "<<mu.pt()<<std::endl;
	//return 0;
	//}
	//}

  //for (const auto &ele : *eles)
	//{
	//std::cout<<"New ELE "<<ele.pt() << " " << ele.electronID("mvaEleID-Fall17-iso-V1-wp80") << " " <<  ele.electronID("mvaEleID-Fall17-noIso-V1-wp80") <<std::endl;
	//if(ele.pt()>60 and ele.electronID("egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wp90"))
	//{
	//std::cout<<"FAIL ele iso pt "<<ele.pt()<<std::endl;
	//return 0;
	//}
	//if(ele.pt()>150 and ele.electronID("mvaEleID-Fall17-noIso-V1-wp90"))
	//{
	//std::cout<<"FAIL ele NONiso pt "<<ele.pt()<<std::endl;
	//return 0;
	//}
	//}
   return 1;

 }




void 
NanoAOD_Filter::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
}

// ------------ method called once each job just before starting event loop  ------------
void 
NanoAOD_Filter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
NanoAOD_Filter::endJob() 
{

}

#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(NanoAOD_Filter);



