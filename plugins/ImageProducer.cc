#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "TLorentzVector.h"

#include <iostream>
#include <fstream>

#include <memory>
#include <iostream>


struct ImageTFCache {
  ImageTFCache() : graphDef(nullptr) {
  }
  //From deepflavour implementation, for consistency 
  std::atomic<tensorflow::GraphDef*> graphDef;
  std::atomic<tensorflow::GraphDef*> graphDefMD;
  //std::atomic<tensorflow::GraphDef*> graphDefPho;
  std::atomic<tensorflow::GraphDef*> graphDefPhoMD;
  std::atomic<tensorflow::GraphDef*> graphDefW;
  std::atomic<tensorflow::GraphDef*> graphDefWMD;
  std::atomic<tensorflow::GraphDef*> graphDefHMD;
  std::atomic<tensorflow::GraphDef*> graphDefHflessMD;
  //std::atomic<tensorflow::GraphDef*> graphDefHfonlyMD;
  std::atomic<tensorflow::GraphDef*> graphDefZMD;
  std::atomic<tensorflow::GraphDef*> graphDefZflessMD;
  //std::atomic<tensorflow::GraphDef*> graphDefZfonlyMD;
  std::atomic<tensorflow::GraphDef*> graphDefWWMD;
  std::atomic<tensorflow::GraphDef*> graphDefWWlepMD;
  std::atomic<tensorflow::GraphDef*> graphDefHWWMD;
  std::atomic<tensorflow::GraphDef*> graphDefHWWlepMD;
  //std::atomic<tensorflow::GraphDef*> graphDefHOT;
  std::atomic<tensorflow::GraphDef*> graphDefMDHOT;
  std::atomic<tensorflow::GraphDef*> graphDefWWMDHOT;
  std::atomic<tensorflow::GraphDef*> graphDefWWlepMDHOT;
  std::atomic<tensorflow::GraphDef*> graphDefHWWMDHOT;
  std::atomic<tensorflow::GraphDef*> graphDefHWWlepMDHOT;

};

class ImageProducer : public edm::stream::EDProducer<edm::GlobalCache<ImageTFCache>> {

  public:
    explicit ImageProducer(const edm::ParameterSet&, const ImageTFCache*);
    ~ImageProducer() override;
    double runtflow(tensorflow::Session* ,tensorflow::Tensor ,tensorflow::Tensor, uint );
    double runtflow(tensorflow::GraphDef* graphDef ,tensorflow::Tensor ,tensorflow::Tensor, uint );
    double principal_axis(const std::vector<std::vector<float>> &);
    static void fillDescriptions(edm::ConfigurationDescriptions&);
    static void globalEndJob(const ImageTFCache*);
    static std::unique_ptr<ImageTFCache> initializeGlobalCache(const edm::ParameterSet&);

  private:
    void beginStream(edm::StreamID) override {}
    void produce(edm::Event&, const edm::EventSetup&) override;
    void endStream() override {}
    ImageTFCache* cache_;
    tensorflow::Session* tfsession_;
    tensorflow::Session* tfsessionMD_;
    //tensorflow::Session* tfsessionPho_;
    tensorflow::Session* tfsessionPhoMD_;
    tensorflow::Session* tfsessionW_;
    tensorflow::Session* tfsessionWMD_;
    tensorflow::Session* tfsessionHMD_;
    tensorflow::Session* tfsessionHflessMD_;
    //tensorflow::Session* tfsessionHfonlyMD_;
    tensorflow::Session* tfsessionZMD_;
    tensorflow::Session* tfsessionZflessMD_;
    //tensorflow::Session* tfsessionZfonlyMD_;
    tensorflow::Session* tfsessionWWMD_;
    tensorflow::Session* tfsessionWWlepMD_;
    tensorflow::Session* tfsessionHWWMD_;
    tensorflow::Session* tfsessionHWWlepMD_;
    //tensorflow::Session* tfsessionHOT_;
    tensorflow::Session* tfsessionMDHOT_;
    tensorflow::Session* tfsessionWWMDHOT_;
    tensorflow::Session* tfsessionWWlepMDHOT_;
    tensorflow::Session* tfsessionHWWMDHOT_;
    tensorflow::Session* tfsessionHWWlepMDHOT_;


    bool isHotVR;

    const edm::EDGetTokenT<edm::View<pat::Jet>> src_;
    const edm::EDGetTokenT<edm::View<pat::Jet>> sj_;
    std::string sdmcoll_;

    edm::FileInPath pb_path_;
    edm::FileInPath pb_pathMD_;
    //edm::FileInPath pb_pathPho_;
    edm::FileInPath pb_pathPhoMD_;
    edm::FileInPath pb_pathW_;
    edm::FileInPath pb_pathWMD_;
    edm::FileInPath pb_pathHMD_;
    edm::FileInPath pb_pathHflessMD_;
    //edm::FileInPath pb_pathHfonlyMD_;
    edm::FileInPath pb_pathZMD_;
    edm::FileInPath pb_pathZflessMD_;
    //edm::FileInPath pb_pathZfonlyMD_;
    edm::FileInPath pb_pathWWMD_;
    edm::FileInPath pb_pathWWlepMD_;
    edm::FileInPath pb_pathHWWMD_;
    edm::FileInPath pb_pathHWWlepMD_;
    //edm::FileInPath pb_pathHOT_;
    edm::FileInPath pb_pathMDHOT_;
    edm::FileInPath pb_pathWWMDHOT_;
    edm::FileInPath pb_pathWWlepMDHOT_;
    edm::FileInPath pb_pathHWWMDHOT_;
    edm::FileInPath pb_pathHWWlepMDHOT_;
    std::string extex_;
    bool isHotVR_;
    uint nsubs;
};

ImageProducer::ImageProducer(const edm::ParameterSet& iConfig,  const ImageTFCache* cache)
: 
   tfsession_(nullptr)
,  tfsessionMD_(nullptr)
//,  tfsessionPho_(nullptr)
,  tfsessionPhoMD_(nullptr)
,  tfsessionW_(nullptr)
,  tfsessionWMD_(nullptr)
,  tfsessionHMD_(nullptr)
,  tfsessionHflessMD_(nullptr)
//,  tfsessionHfonlyMD_(nullptr)
,  tfsessionZMD_(nullptr)
,  tfsessionZflessMD_(nullptr)
//,  tfsessionZfonlyMD_(nullptr)
,  tfsessionWWMD_(nullptr)
,  tfsessionWWlepMD_(nullptr)
,  tfsessionHWWMD_(nullptr)
,  tfsessionHWWlepMD_(nullptr)
//,  tfsessionHOT_(nullptr)
,  tfsessionMDHOT_(nullptr)
,  tfsessionWWMDHOT_(nullptr)
,  tfsessionWWlepMDHOT_(nullptr)
,  tfsessionHWWMDHOT_(nullptr)
,  tfsessionHWWlepMDHOT_(nullptr)
,  src_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("src")))
,  sj_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("sj")))
,  sdmcoll_(iConfig.getParameter<std::string>("sdmcoll"))
,  extex_(iConfig.getParameter<std::string>("extex"))
,  isHotVR_(iConfig.getParameter<bool>("isHotVR"))
{
  
  produces<pat::JetCollection>();

  tensorflow::SessionOptions sessionOptions;
  const ImageTFCache* cache_=cache;

  tfsession_ = tensorflow::createSession(cache_->graphDef,sessionOptions);
  tfsessionMD_ = tensorflow::createSession(cache_->graphDefMD,sessionOptions);
  //tfsessionPho_ = tensorflow::createSession(cache_->graphDefPho,sessionOptions);
  tfsessionPhoMD_ = tensorflow::createSession(cache_->graphDefPhoMD,sessionOptions);
  tfsessionW_ = tensorflow::createSession(cache_->graphDefW,sessionOptions);
  tfsessionWMD_ = tensorflow::createSession(cache_->graphDefWMD,sessionOptions);
  tfsessionHMD_ = tensorflow::createSession(cache_->graphDefHMD,sessionOptions);
  tfsessionHflessMD_ = tensorflow::createSession(cache_->graphDefHflessMD,sessionOptions);
  //tfsessionHfonlyMD_ = tensorflow::createSession(cache_->graphDefHfonlyMD,sessionOptions);
  tfsessionZMD_ = tensorflow::createSession(cache_->graphDefZMD,sessionOptions);
  tfsessionZflessMD_ = tensorflow::createSession(cache_->graphDefZflessMD,sessionOptions);
  //tfsessionZfonlyMD_ = tensorflow::createSession(cache_->graphDefZfonlyMD,sessionOptions);
  tfsessionWWMD_ = tensorflow::createSession(cache_->graphDefWWMD,sessionOptions);
  tfsessionWWlepMD_ = tensorflow::createSession(cache_->graphDefWWlepMD,sessionOptions);
  tfsessionHWWMD_ = tensorflow::createSession(cache_->graphDefHWWMD,sessionOptions);
  tfsessionHWWlepMD_ = tensorflow::createSession(cache_->graphDefHWWlepMD,sessionOptions);
  //tfsessionHOT_ = tensorflow::createSession(cache_->graphDefHOT,sessionOptions);
  tfsessionMDHOT_ = tensorflow::createSession(cache_->graphDefMDHOT,sessionOptions);
  tfsessionWWMDHOT_ = tensorflow::createSession(cache_->graphDefWWMDHOT,sessionOptions);
  tfsessionWWlepMDHOT_ = tensorflow::createSession(cache_->graphDefWWlepMDHOT,sessionOptions);
  tfsessionHWWMDHOT_ = tensorflow::createSession(cache_->graphDefHWWMDHOT,sessionOptions);
  tfsessionHWWlepMDHOT_ = tensorflow::createSession(cache_->graphDefHWWlepMDHOT,sessionOptions);
  isHotVR=isHotVR_;


  //if(iConfig.getParameter<edm::InputTag>("src").label()=="slimmedJetsAK8")sdmcoll="ak8PFJetsCHSSoftDropMass";

}
ImageProducer::~ImageProducer(){
  if (tfsession_ != nullptr) tensorflow::closeSession(tfsession_);
  if (tfsessionMD_ != nullptr)  tensorflow::closeSession(tfsessionMD_);
  //if (tfsessionPho_ != nullptr)  tensorflow::closeSession(tfsessionPho_);
  if (tfsessionPhoMD_ != nullptr)  tensorflow::closeSession(tfsessionPhoMD_);
  if (tfsessionW_ != nullptr)  tensorflow::closeSession(tfsessionW_);
  if (tfsessionWMD_ != nullptr)  tensorflow::closeSession(tfsessionWMD_);
  if (tfsessionHMD_ != nullptr)  tensorflow::closeSession(tfsessionHMD_);
  if (tfsessionHflessMD_ != nullptr)  tensorflow::closeSession(tfsessionHflessMD_);
  //if (tfsessionHfonlyMD_ != nullptr)  tensorflow::closeSession(tfsessionHfonlyMD_);
  if (tfsessionZMD_ != nullptr)  tensorflow::closeSession(tfsessionZMD_);
  if (tfsessionZflessMD_ != nullptr)  tensorflow::closeSession(tfsessionZflessMD_);
  //if (tfsessionZfonlyMD_ != nullptr)  tensorflow::closeSession(tfsessionZfonlyMD_);
  if (tfsessionWWMD_ != nullptr)  tensorflow::closeSession(tfsessionWWMD_);
  if (tfsessionWWlepMD_ != nullptr)  tensorflow::closeSession(tfsessionWWlepMD_);
  if (tfsessionHWWMD_ != nullptr)  tensorflow::closeSession(tfsessionHWWMD_);
  if (tfsessionHWWlepMD_ != nullptr)  tensorflow::closeSession(tfsessionHWWlepMD_);
  //if (tfsessionHOT_ != nullptr)  tensorflow::closeSession(tfsessionHOT_);
  if (tfsessionMDHOT_ != nullptr)  tensorflow::closeSession(tfsessionMDHOT_);
  if (tfsessionWWMDHOT_ != nullptr)  tensorflow::closeSession(tfsessionWWMDHOT_);
  if (tfsessionWWlepMDHOT_ != nullptr)  tensorflow::closeSession(tfsessionWWlepMDHOT_);
  if (tfsessionHWWMDHOT_ != nullptr)  tensorflow::closeSession(tfsessionHWWMDHOT_);
  if (tfsessionHWWlepMDHOT_ != nullptr)  tensorflow::closeSession(tfsessionHWWlepMDHOT_);
}

void ImageProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


std::unique_ptr<ImageTFCache> ImageProducer::initializeGlobalCache(
  const edm::ParameterSet& iConfig)
{
  tensorflow::setLogging("3");
  ImageTFCache* cache = new ImageTFCache();

  cache->graphDef = tensorflow::loadGraphDef(iConfig.getUntrackedParameter<edm::FileInPath>("pb_path").fullPath());
  cache->graphDefMD = tensorflow::loadGraphDef(iConfig.getUntrackedParameter<edm::FileInPath>("pb_pathMD").fullPath());
  //cache->graphDefPho = tensorflow::loadGraphDef(iConfig.getUntrackedParameter<edm::FileInPath>("pb_pathPho").fullPath());
  cache->graphDefPhoMD = tensorflow::loadGraphDef(iConfig.getUntrackedParameter<edm::FileInPath>("pb_pathPhoMD").fullPath());
  cache->graphDefW = tensorflow::loadGraphDef(iConfig.getUntrackedParameter<edm::FileInPath>("pb_pathW").fullPath());
  cache->graphDefWMD = tensorflow::loadGraphDef(iConfig.getUntrackedParameter<edm::FileInPath>("pb_pathWMD").fullPath());
  cache->graphDefHMD = tensorflow::loadGraphDef(iConfig.getUntrackedParameter<edm::FileInPath>("pb_pathHMD").fullPath());
  cache->graphDefHflessMD = tensorflow::loadGraphDef(iConfig.getUntrackedParameter<edm::FileInPath>("pb_pathHflessMD").fullPath());
  //cache->graphDefHfonlyMD = tensorflow::loadGraphDef(iConfig.getUntrackedParameter<edm::FileInPath>("pb_pathHfonlyMD").fullPath());
  cache->graphDefZMD = tensorflow::loadGraphDef(iConfig.getUntrackedParameter<edm::FileInPath>("pb_pathZMD").fullPath());
  cache->graphDefZflessMD = tensorflow::loadGraphDef(iConfig.getUntrackedParameter<edm::FileInPath>("pb_pathZflessMD").fullPath());
  //cache->graphDefZfonlyMD = tensorflow::loadGraphDef(iConfig.getUntrackedParameter<edm::FileInPath>("pb_pathZfonlyMD").fullPath());
  cache->graphDefWWMD = tensorflow::loadGraphDef(iConfig.getUntrackedParameter<edm::FileInPath>("pb_pathWWMD").fullPath());
  cache->graphDefWWlepMD = tensorflow::loadGraphDef(iConfig.getUntrackedParameter<edm::FileInPath>("pb_pathWWlepMD").fullPath());
  cache->graphDefHWWMD = tensorflow::loadGraphDef(iConfig.getUntrackedParameter<edm::FileInPath>("pb_pathHWWMD").fullPath());
  cache->graphDefHWWlepMD = tensorflow::loadGraphDef(iConfig.getUntrackedParameter<edm::FileInPath>("pb_pathHWWlepMD").fullPath());
  //cache->graphDefHOT = tensorflow::loadGraphDef(iConfig.getUntrackedParameter<edm::FileInPath>("pb_pathHOT").fullPath());
  
  cache->graphDefMDHOT = tensorflow::loadGraphDef(iConfig.getUntrackedParameter<edm::FileInPath>("pb_pathMDHOT").fullPath());
  cache->graphDefWWMDHOT = tensorflow::loadGraphDef(iConfig.getUntrackedParameter<edm::FileInPath>("pb_pathWWMDHOT").fullPath());
  cache->graphDefWWlepMDHOT = tensorflow::loadGraphDef(iConfig.getUntrackedParameter<edm::FileInPath>("pb_pathWWlepMDHOT").fullPath());
  cache->graphDefHWWMDHOT = tensorflow::loadGraphDef(iConfig.getUntrackedParameter<edm::FileInPath>("pb_pathHWWMDHOT").fullPath());
  cache->graphDefHWWlepMDHOT = tensorflow::loadGraphDef(iConfig.getUntrackedParameter<edm::FileInPath>("pb_pathHWWlepMDHOT").fullPath());

  return std::unique_ptr<ImageTFCache>(cache);
}


double ImageProducer::runtflow(tensorflow::Session* tfsession_,tensorflow::Tensor input_image,tensorflow::Tensor input_b,uint NQCDs)
	{

	//tensorflow::Status status1 = tfsession_->Reset({ { "input_1", input_image }, {"input_2", input_b} },{ "k2tfout_0" }, {}, &tfoutput);
	std::vector<tensorflow::Tensor> tfoutput;
  	tensorflow::Status status = tfsession_->Run({ { "input_1", input_image }, {"input_2", input_b} },{ "k2tfout_0" }, {}, &tfoutput);

    	//tensorflow::run(tfsession_, { { "input_1", input_image }, {"input_2", input_b} },{ "k2tfout_0" }, {}, &tfoutput);
	if (!status.ok()) 
		{
		std::cout << "Tensorflow Failed:" << std::endl;
  		std::cout << status.ToString() << "\n";
  		return -1.0;
		}	
        float result_top = 0.0;
        float result_qcd = 0.0;
	
	for(uint i=0;i <  tfoutput[0].flat<float>().size();++i) 
		{
		if (tfoutput[0].flat<float>()(i)==1.0)
			{
			//std::cout << "ERR index:" << i<<std::endl;
			return 0.0;
			}
		//std::cout<<i<<": "<<tfoutput[0].flat<float>()(i)<<std::endl;
		if(i<NQCDs) result_qcd += tfoutput[0].flat<float>()(i);
		else result_top += tfoutput[0].flat<float>()(i);
		}
	if ((result_top+result_qcd)==0.0)
		{
		std::cout << "ERR result_top:" << result_top<<" result_qcd:"<<result_qcd<<std::endl;
		return 0.0;
		}
	return result_top/(result_top+result_qcd);
	 
	}


double ImageProducer::principal_axis(const std::vector<std::vector<float>> & partlist)
	{
  	double tan_theta=0.0;
	double M11=0.0;
	double M20=0.0;
	double M02=0.0;
	for(uint i=0;i < partlist[0].size();++i)
		{
		M11 += partlist[0][i]*partlist[1][i]*partlist[2][i];
		M20 += partlist[0][i]*partlist[1][i]*partlist[1][i];
		M02 += partlist[0][i]*partlist[2][i]*partlist[2][i];
		}
  	double denom=(M20-M02+std::sqrt(4*M11*M11+(M20-M02)*(M20-M02)));
  	if(denom!=0.0) tan_theta=2.0*M11/denom;
  	return tan_theta;
	}


void ImageProducer::globalEndJob(const ImageTFCache* cache)
{
  delete cache;  
}

void ImageProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<edm::View<pat::Jet>> jets;
  edm::Handle<edm::View<pat::Jet>> subjets;

  iEvent.getByToken(src_, jets);
  iEvent.getByToken(sj_, subjets);

  TLorentzVector curtlv;
  TLorentzVector sublv;

  int jindex=0;

  std::vector<float> itopdisc = {};
  std::vector<float> itopdiscMD = {};
  //std::vector<float> itopdiscPho = {};
  std::vector<float> itopdiscPhoMD = {};
  std::vector<float> itopdiscW = {};
  std::vector<float> itopdiscWMD = {};
  std::vector<float> itopdiscHMD = {};
  std::vector<float> itopdiscHflessMD = {};
  //std::vector<float> itopdiscHfonlyMD = {};
  std::vector<float> itopdiscZMD = {};
  std::vector<float> itopdiscZflessMD = {};
  //std::vector<float> itopdiscZfonlyMD = {};
  std::vector<float> itopdiscWWMD = {};
  std::vector<float> itopdiscWWlepMD = {};
  std::vector<float> itopdiscHWWMD = {};
  std::vector<float> itopdiscHWWlepMD = {};
  //std::vector<float> itopdiscHOT = {};
  std::vector<float> itopdiscMDHOT = {};
  std::vector<float> itopdiscWWMDHOT = {};
  std::vector<float> itopdiscWWlepMDHOT = {};
  std::vector<float> itopdiscHWWMDHOT = {};
  std::vector<float> itopdiscHWWlepMDHOT = {};
  std::vector<float> gmasses = {};

  if(isHotVR) nsubs=4;
  else nsubs=2;
  int ntopinit = -1;
  for (const auto &AK8pfjet : *jets)
	{
	ntopinit+=1;


    	if(isHotVR)
		{
	  	//itopdiscHOT.push_back(-10.0);
	  	itopdiscMDHOT.push_back(-10.0);
	  	itopdiscWWMDHOT.push_back(-10.0);
	  	itopdiscWWlepMDHOT.push_back(-10.0);
	  	itopdiscHWWMDHOT.push_back(-10.0);
	  	itopdiscHWWlepMDHOT.push_back(-10.0);



		}
	else 
		{  	
		itopdisc.push_back(-10.0);
	  	itopdiscMD.push_back(-10.0);
	  	//itopdiscPho.push_back(-10.0);
	  	itopdiscPhoMD.push_back(-10.0);
	  	itopdiscW.push_back(-10.0);
	  	itopdiscWMD.push_back(-10.0);
	  	itopdiscHMD.push_back(-10.0);
	  	itopdiscHflessMD.push_back(-10.0);
	  	//itopdiscHfonlyMD.push_back(-10.0);
	  	itopdiscZMD.push_back(-10.0);
	  	itopdiscZflessMD.push_back(-10.0);
	  	//itopdiscZfonlyMD.push_back(-10.0);
	  	itopdiscWWMD.push_back(-10.0);
	  	itopdiscWWlepMD.push_back(-10.0);
	  	itopdiscHWWMD.push_back(-10.0);
	  	itopdiscHWWlepMD.push_back(-10.0);
		}

  	gmasses.push_back(-10.0);

        TLorentzVector curtlv;
	curtlv.SetPtEtaPhiM(AK8pfjet.pt(),AK8pfjet.eta(),AK8pfjet.phi(),AK8pfjet.mass());


	float mergeval=0.8;
	if(isHotVR) mergeval=1.2;

	int ndau = AK8pfjet.numberOfDaughters();
        std::vector<pat::PackedCandidate> allpf;

	std::vector<std::vector<float>> partlist = {{},{},{},{},{},{},{},{}};
	std::vector<float> sjlist = {};
	
	double fullint = 0.0;
        double etac=0;
        double phic=0;

	int idaufill = 0;

	for(int idau=0;idau<ndau;idau++)
		{
	        const pat::PackedCandidate* lPack = dynamic_cast<const pat::PackedCandidate *>(AK8pfjet.daughter(idau) );

	  	if (lPack != nullptr)
			{
			float dphi = reco::deltaPhi(lPack->phi(),curtlv.Phi());

        		TLorentzVector pfclv;
			pfclv.SetPtEtaPhiM(lPack->pt(),lPack->eta(),curtlv.Phi()+dphi,lPack->mass());
			if ((pfclv.Pt()<=1.0) and (lPack->charge()!=0)) continue;
			double funcval =(6.62733e-02) + (2.63732e+02)*(1.0/curtlv.Perp());
			double drcorval = 0.6/(funcval);

			int pfcflav = std::abs(lPack->pdgId());

			double neweta = pfclv.Eta()+(pfclv.Eta()-curtlv.Eta())*(drcorval-1.0);
			double newphi = pfclv.Phi()+(dphi)*(drcorval-1.0);

			double newdetafj = (pfclv.Eta()-curtlv.Eta())*drcorval;
			double newdphifj = (dphi)*drcorval;

			if(std::sqrt(newdphifj*newdphifj+newdetafj*newdetafj)>1.6) continue;

			partlist[0].push_back(lPack->pt());
			partlist[1].push_back(neweta);
			partlist[2].push_back(newphi);

			fullint+=partlist[0][idaufill];
			etac += partlist[0][idaufill]*partlist[1][idaufill];
			phic += partlist[0][idaufill]*partlist[2][idaufill];

			if(pfcflav==13)partlist[3].push_back(lPack->pt());
			else partlist[3].push_back(0.0);
			if(pfcflav==11)partlist[4].push_back(lPack->pt());
			else partlist[4].push_back(0.0);
			if(pfcflav==211)partlist[5].push_back(lPack->pt());
			else partlist[5].push_back(0.0);
			if(pfcflav==22)partlist[6].push_back(lPack->pt());
			else partlist[6].push_back(0.0);
			if(pfcflav==130)partlist[7].push_back(lPack->pt());
			else partlist[7].push_back(0.0);
			idaufill+=1;
			}
		}
	TLorentzVector gjet;

        std::vector<pat::Jet> sjvec;
        std::vector<pat::Jet> sjvecmatch;
	std::vector<int>lepblacks{5,11};
  	for (const auto &subjet : *subjets)
		{
	       

		sjvec.push_back(subjet);
		sublv.SetPtEtaPhiM(subjet.pt(),subjet.eta(),subjet.phi(),subjet.mass());
		
		if (sublv.DeltaR(curtlv)>mergeval || sjlist.size()>=(nsubs*6)) continue;
		sjvecmatch.push_back(subjet);
		if(sjlist.size()==0)gjet=sublv;
		else gjet+=sublv;
		sjlist.push_back(subjet.bDiscriminator("pfDeepFlavourJetTags:probb"));
		sjlist.push_back(subjet.bDiscriminator("pfDeepFlavourJetTags:probbb"));
		sjlist.push_back(subjet.bDiscriminator("pfDeepFlavourJetTags:probuds"));
		sjlist.push_back(subjet.bDiscriminator("pfDeepFlavourJetTags:probg"));
		sjlist.push_back(subjet.bDiscriminator("pfDeepFlavourJetTags:probc"));
		sjlist.push_back(subjet.bDiscriminator("pfDeepFlavourJetTags:problepb"));
		}
	float gmass = 0.0;

	if(isHotVR)	
		{
		//nsjh = sjvecmatch.size();
		//mmin = 0.;
		//fpt = 0.;
		//if (nsjh >= 3)	{
		//	  	fpt = sjvecmatch[0].pt() / AK8pfjet.pt();
		//	      	mmin = std::min({
		//		(sjvecmatch[0].p4()+sjvecmatch[1].p4()).mass(),
		//		(sjvecmatch[0].p4()+sjvecmatch[2].p4()).mass(),
		//		(sjvecmatch[1].p4()+sjvecmatch[2].p4()).mass(),
		//	      	});
		//  		} 
		if (sjvecmatch.size()>0) gmass=gjet.M();
		}
	else	gmass=std::max(float(0.0),AK8pfjet.userFloat(sdmcoll_));


	uint sjlsize=sjlist.size();

	for(uint isj=0;isj<(nsubs*6-sjlsize);isj++) sjlist.push_back(0.);
	uint DBindex=sjlist.size();

	sjlist.push_back(AK8pfjet.bDiscriminator("pfMassIndependentDeepDoubleBvLJetTags:probHbb"));

        sjlist.push_back(gmass/172.0);
        sjlist.push_back(AK8pfjet.pt()/2000.0);
	gmasses[ntopinit]=gmass;

	


        int npoints=38;
	std::vector<int> ietalist = {};
	std::vector<int> iphilist = {};



	//centering and rotating 
	etac/=fullint;
	phic/=fullint;
	for(uint i=0;i < partlist[0].size();++i)
		{
	        partlist[1][i] -= etac;
	        partlist[2][i] -= phic;
		}

	double tan_theta=principal_axis(partlist); 
	double DReta=1.6;
	double DRphi=1.6;
	for(uint i=0;i < partlist[0].size();++i)
		{
		double Reta = partlist[1][i]*std::cos(std::atan(tan_theta))+partlist[2][i]*std::sin(std::atan(tan_theta));
		double Rphi = -1.0*partlist[1][i]*std::sin(std::atan(tan_theta))+partlist[2][i]*std::cos(std::atan(tan_theta));

		partlist[1][i] = Reta;
		partlist[2][i] = Rphi;

		ietalist.push_back(int((partlist[1][i]+DReta)/(2.0*DReta/float(npoints-1))));
		iphilist.push_back(int((partlist[2][i]+DRphi)/(2.0*DRphi/float(npoints-1))));
		}

  	uint ncolors=6;
  
	std::vector<std::vector<std::vector<double>>> grid(37,std::vector<std::vector<double>>(37,std::vector<double>(ncolors,0.0)));
	std::vector<std::pair<std::vector<uint>,std::vector<double>>> indexedimage;


	//normalization and digitization
	for(uint i=0;i < partlist[0].size();++i)
		{
	        if((ietalist[i]>=37) or (iphilist[i]>=37) or (ietalist[i]<=0) or (iphilist[i]<=0))continue;
		int filldex=0;

		for(uint j=0;j < partlist.size();++j)
			{
			if(((j>2) or (j==0))) //1 and 2 are eta,phi
				{
				grid[ietalist[i]][iphilist[i]][filldex]+=partlist[j][i]/fullint;
				filldex+=1;
				}
			}				
		}

	for(uint i=0;i < grid.size();++i)
		{
		for(uint j=0;j < grid[i].size();++j)
			{
				if(grid[i][j][0]>0.0000000001)
					{
					std::pair<std::vector<uint>,std::vector<double>> elem;
					elem.first = {i,j};
					for(uint k=0;k < grid[i][j].size();++k)elem.second.push_back(grid[i][j][k]);					
					indexedimage.push_back(elem);
					}
			}
		}
	
	//flipping horiz and vert
	uint half_img=(npoints-2)/2;
	float left_sum=0.0;
	float right_sum=0.0;
	for(uint i=0;i < indexedimage.size();++i)
		{
		if(indexedimage[i].first[0]<half_img)left_sum+=indexedimage[i].second[0];
		if(indexedimage[i].first[0]>half_img)right_sum+=indexedimage[i].second[0];
		}
	if(left_sum<right_sum)
		{
		for(uint i=0;i < indexedimage.size();++i)indexedimage[i].first={npoints-2-indexedimage[i].first[0],indexedimage[i].first[1]};	
		}

	float lower_sum=0.0;
	float upper_sum=0.0;
	for(uint i=0;i < indexedimage.size();++i)
		{
		if(indexedimage[i].first[1]>half_img)lower_sum+=indexedimage[i].second[0];
		if(indexedimage[i].first[1]<half_img)upper_sum+=indexedimage[i].second[0];
		}

	if(lower_sum<upper_sum)
		{
		for(uint i=0;i < indexedimage.size();++i)indexedimage[i].first={indexedimage[i].first[0],npoints-2-indexedimage[i].first[1]};		
		}

	//convert scalars to tensorflow
	tensorflow::Tensor input_b(tensorflow::DT_FLOAT, { 1, int(sjlist.size()) });
	float* d = input_b.flat<float>().data();

	tensorflow::Tensor input_nodoubleb(tensorflow::DT_FLOAT, { 1, int(sjlist.size())-1 });
	float* dndb = input_nodoubleb.flat<float>().data();

	tensorflow::Tensor input_nolep(tensorflow::DT_FLOAT, { 1, int(sjlist.size() - lepblacks.size()) });
	float* dnl = input_nolep.flat<float>().data();
	
	tensorflow::Tensor input_fless(tensorflow::DT_FLOAT, { 1, 2 });
	float* dfless = input_fless.flat<float>().data();


	for(uint i=0;i < sjlist.size();++i,++d)
		{
		
		*d = sjlist[i];	

		if (i!=DBindex) 
			{
			*dndb = sjlist[i];
			++dndb;
			}
		if (!(std::find(lepblacks.begin(), lepblacks.end(), i) != lepblacks.end())) 
			{
			*dnl = sjlist[i];
			++dnl;
			}
		if((sjlist.size()-i)<3) 
			{
			*dfless = sjlist[i];
			++dfless;
			}
		}	
	//convert image to tensorflow.  first create tensor of zeros, then fill.  Probably not optimal quite yet
	tensorflow::Tensor input_image(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1,37, 37 , ncolors }));
	auto input_map = input_image.tensor<float, 4>();

	tensorflow::Tensor input_image_nolep(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1,37, 37 , 3 }));
	auto input_map_nolep = input_image_nolep.tensor<float, 4>();

	for(uint i=0;i < 37;++i)
		{
		for(uint j=0;j < 37;++j)
			{
			for(uint k=0;k < ncolors;++k)
				{
				input_map(0,i,j,k) = 0.0;
				if(k>2)input_map_nolep(0,i,j,k-3) = 0.0;
				}
			}
			
		}
	
	for(uint i=0;i < indexedimage.size();++i)
		{
		for(uint j=0;j < indexedimage[i].second.size();++j)
			{
				input_map(0,indexedimage[i].first[0], indexedimage[i].first[1], j) = indexedimage[i].second[j];
				if(j>2)input_map_nolep(0,indexedimage[i].first[0], indexedimage[i].first[1], j-3) = indexedimage[i].second[j];
			}	
		}



	
	//Actually run tensorflow

    	if(isHotVR)
		{
		itopdiscMDHOT[jindex]=runtflow(tfsessionMDHOT_,input_image,input_nodoubleb,1);
		itopdiscWWMDHOT[jindex]=runtflow(tfsessionWWMDHOT_,input_image,input_nodoubleb,1);
		itopdiscWWlepMDHOT[jindex]=runtflow(tfsessionWWlepMDHOT_,input_image,input_nodoubleb,1);
		itopdiscHWWMDHOT[jindex]=runtflow(tfsessionHWWMDHOT_,input_image,input_nodoubleb,1);
		itopdiscHWWlepMDHOT[jindex]=runtflow(tfsessionHWWlepMDHOT_,input_image,input_nodoubleb,1);
		}
    	else
		{        
		itopdisc[jindex]=runtflow(tfsession_,input_image,input_nodoubleb,4);
		itopdiscMD[jindex]=runtflow(tfsessionMD_,input_image,input_nodoubleb,4);
		itopdiscPhoMD[jindex]=runtflow(tfsessionPhoMD_,input_image_nolep,input_nolep,4);
		itopdiscW[jindex]=runtflow(tfsessionW_,input_image,input_nodoubleb,4);
		itopdiscWMD[jindex]=runtflow(tfsessionWMD_,input_image,input_nodoubleb,4);
		itopdiscHflessMD[jindex]=runtflow(tfsessionHflessMD_,input_image,input_fless,4);
		itopdiscHMD[jindex]=runtflow(tfsessionHMD_,input_image_nolep,input_nolep,4);
		//itopdiscHfonlyMD[jindex]=runtflow(tfsessionHfonlyMD_,input_image,input_nolep,4);
		itopdiscZflessMD[jindex]=runtflow(tfsessionZflessMD_,input_image,input_fless,4);
		itopdiscZMD[jindex]=runtflow(tfsessionZMD_,input_image_nolep,input_nolep,4);
		//itopdiscZfonlyMD[jindex]=runtflow(tfsessionZfonlyMD_,input_image,input_nolep,4);
		itopdiscWWMD[jindex]=runtflow(tfsessionWWMD_,input_image,input_nodoubleb,4);
		itopdiscWWlepMD[jindex]=runtflow(tfsessionWWlepMD_,input_image,input_nodoubleb,4);
		itopdiscHWWMD[jindex]=runtflow(tfsessionHWWMD_,input_image,input_nodoubleb,4);
		itopdiscHWWlepMD[jindex]=runtflow(tfsessionHWWlepMD_,input_image,input_nodoubleb,4);
		}
	jindex+=1;

       
  }


  //Add to jet userfloats
  auto outputs = std::make_unique<pat::JetCollection>();
  jindex=0;

  for (const auto &jet : *jets){
    pat::Jet newJet(jet);
    if(isHotVR)
		{
    		//newJet.addUserFloat("Image"+extex_+":top", itopdisc[jindex]);
    		newJet.addUserFloat("ImageMD"+extex_+":top", itopdiscMDHOT[jindex]);
    		newJet.addUserFloat("ImageMD"+extex_+":ww", itopdiscWWMDHOT[jindex]);
    		newJet.addUserFloat("ImageMD"+extex_+":wwlep", itopdiscWWlepMDHOT[jindex]);
    		newJet.addUserFloat("ImageMD"+extex_+":hww", itopdiscHWWMDHOT[jindex]);
    		newJet.addUserFloat("ImageMD"+extex_+":hwwlep", itopdiscHWWlepMDHOT[jindex]);

		}
    else
		{
	    	newJet.addUserFloat("Image"+extex_+":top", itopdisc[jindex]);
	    	newJet.addUserFloat("ImageMD"+extex_+":top", itopdiscMD[jindex]);
	    	//newJet.addUserFloat("Image"+extex_+":pho", itopdiscPho[jindex]);
	    	newJet.addUserFloat("ImageMD"+extex_+":pho", itopdiscPhoMD[jindex]);
	    	newJet.addUserFloat("Image"+extex_+":w", itopdiscW[jindex]);
	    	newJet.addUserFloat("ImageMD"+extex_+":w", itopdiscWMD[jindex]);
	    	newJet.addUserFloat("ImageMD"+extex_+":h", itopdiscHMD[jindex]);
	    	newJet.addUserFloat("ImageMD"+extex_+":hfless", itopdiscHflessMD[jindex]);
	    	//newJet.addUserFloat("ImageMD"+extex_+":hfonly", itopdiscHfonlyMD[jindex]);
	    	newJet.addUserFloat("ImageMD"+extex_+":z", itopdiscZMD[jindex]);
	    	newJet.addUserFloat("ImageMD"+extex_+":zfless", itopdiscZflessMD[jindex]);
	    	//newJet.addUserFloat("ImageMD"+extex_+":zfonly", itopdiscZfonlyMD[jindex]);
	    	newJet.addUserFloat("ImageMD"+extex_+":ww", itopdiscWWMD[jindex]);
	    	newJet.addUserFloat("ImageMD"+extex_+":wwlep", itopdiscWWlepMD[jindex]);
	    	newJet.addUserFloat("ImageMD"+extex_+":hww", itopdiscHWWMD[jindex]);
	    	newJet.addUserFloat("ImageMD"+extex_+":hwwlep", itopdiscHWWlepMD[jindex]);
		}

    newJet.addUserFloat("Image"+extex_+":mass", gmasses[jindex]);

    outputs->push_back(newJet);
    jindex+=1;

  }
  // put into the event
  iEvent.put(std::move(outputs));
}

//define this as a plug-in
DEFINE_FWK_MODULE(ImageProducer);
