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
};

class ImageProducer : public edm::stream::EDProducer<edm::GlobalCache<ImageTFCache>> {

  public:
    explicit ImageProducer(const edm::ParameterSet&, const ImageTFCache*);
    ~ImageProducer() override;
    double principal_axis(const std::vector<std::vector<float>> &);
    static void fillDescriptions(edm::ConfigurationDescriptions&);
    static void globalEndJob(const ImageTFCache*);
    static std::unique_ptr<ImageTFCache> initializeGlobalCache(const edm::ParameterSet&);

  private:
    void beginStream(edm::StreamID) override {}
    void produce(edm::Event&, const edm::EventSetup&) override;
    void endStream() override {}
    tensorflow::Session* tfsession_;

    const edm::EDGetTokenT<edm::View<pat::Jet>> src_;
    const edm::EDGetTokenT<edm::View<pat::Jet>> sj_;
    std::string sdmcoll_;
    edm::FileInPath pb_path_;
    std::string extex_;

};

ImageProducer::ImageProducer(const edm::ParameterSet& iConfig,  const ImageTFCache* cache)
: 
 tfsession_(nullptr)
, src_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("src")))
, sj_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("sj")))
, sdmcoll_(iConfig.getParameter<std::string>("sdmcoll"))
, extex_(iConfig.getParameter<std::string>("extex"))
{
  
  produces<pat::JetCollection>();

  tensorflow::SessionOptions sessionOptions;
  tfsession_ = tensorflow::createSession(cache->graphDef,sessionOptions);
  //if(iConfig.getParameter<edm::InputTag>("src").label()=="slimmedJetsAK8")sdmcoll="ak8PFJetsCHSSoftDropMass";

}
ImageProducer::~ImageProducer(){
 tensorflow::closeSession(tfsession_);
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
  return std::unique_ptr<ImageTFCache>(cache);
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
  if (cache->graphDef != nullptr) {
    delete cache->graphDef;
  }
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
  int ntopinit = -1;
  //std::cout<<std::endl;
  for (const auto &AK8pfjet : *jets)
	{
	ntopinit+=1;
  	itopdisc.push_back(-10.0);

        TLorentzVector curtlv;
	curtlv.SetPtEtaPhiM(AK8pfjet.pt(),AK8pfjet.eta(),AK8pfjet.phi(),AK8pfjet.mass());

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
        		TLorentzVector pfclv;
			pfclv.SetPtEtaPhiM(lPack->pt(),lPack->eta(),lPack->phi(),lPack->mass());
			if ((pfclv.Pt()<=1.0) and (lPack->charge()!=0)) continue;
			double funcval =(6.62733e-02) + (2.63732e+02)*(1.0/curtlv.Perp());
			double drcorval = 0.6/(funcval);

			int pfcflav = std::abs(lPack->pdgId());

			double neweta = pfclv.Eta()+(pfclv.Eta()-curtlv.Eta())*(drcorval-1.0);
			double newphi = pfclv.Phi()+(pfclv.Phi()-curtlv.Phi())*(drcorval-1.0);

			double newdetafj = (pfclv.Eta()-curtlv.Eta())*drcorval;
			double newdphifj = (pfclv.Phi()-curtlv.Phi())*drcorval;

			if(std::sqrt(newdphifj*newdphifj+newdetafj*newdetafj)>1.6) continue;
			float pw = lPack->puppiWeight();

			partlist[0].push_back(lPack->pt()*pw);
			partlist[1].push_back(neweta);
			partlist[2].push_back(newphi);

			fullint+=partlist[0][idaufill];
			etac += partlist[0][idaufill]*partlist[1][idaufill];
			phic += partlist[0][idaufill]*partlist[2][idaufill];

			if(pfcflav==13)partlist[3].push_back(lPack->pt()*pw);
			else partlist[3].push_back(0.0);
			if(pfcflav==11)partlist[4].push_back(lPack->pt()*pw);
			else partlist[4].push_back(0.0);
			if(pfcflav==211)partlist[5].push_back(lPack->pt()*pw);
			else partlist[5].push_back(0.0);
			if(pfcflav==22)partlist[6].push_back(lPack->pt()*pw);
			else partlist[6].push_back(0.0);
			if(pfcflav==130)partlist[7].push_back(lPack->pt()*pw);
			else partlist[7].push_back(0.0);
			idaufill+=1;
			}
		}

  	for (const auto &subjet : *subjets)
		{
	       
		sublv.SetPtEtaPhiM(subjet.pt(),subjet.eta(),subjet.phi(),subjet.mass());
		if (sublv.DeltaR(curtlv)>0.8 || sjlist.size()>=12) continue;

		sjlist.push_back(subjet.bDiscriminator("pfDeepFlavourJetTags:probb"));
		sjlist.push_back(subjet.bDiscriminator("pfDeepFlavourJetTags:probbb"));
		sjlist.push_back(subjet.bDiscriminator("pfDeepFlavourJetTags:probuds"));
		sjlist.push_back(subjet.bDiscriminator("pfDeepFlavourJetTags:probg"));
		sjlist.push_back(subjet.bDiscriminator("pfDeepFlavourJetTags:probc"));
		sjlist.push_back(subjet.bDiscriminator("pfDeepFlavourJetTags:problepb"));
		}


	uint sjlsize=sjlist.size();
	for(uint isj=0;isj<(12-sjlsize);isj++) sjlist.push_back(0.);
        sjlist.push_back(fabs(AK8pfjet.userFloat(sdmcoll_))/172.0);
        //sjlist.push_back(AK8pfjet.pt());
        //sjlist.push_back(AK8pfjet.eta());
	


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
	std::vector<tensorflow::Tensor> outputs1;
	tensorflow::Tensor input_b(tensorflow::DT_FLOAT, { 1, int(sjlist.size()) });
	float* d = input_b.flat<float>().data();
	uint index=-1;
	for(uint i=0;i < sjlist.size();++i,++d)
		{
		if (i<12)index = (i%6)*2 + i/6;
		else index=i;
		*d = sjlist[index];	
		}	
	
	//convert image to tensorflow.  first create tensor of zeros, then fill.  Probably not optimal quite yet
	tensorflow::Tensor input_image(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1,37, 37 , ncolors }));
	auto input_map = input_image.tensor<float, 4>();

	for(uint i=0;i < 37;++i)
		{
		for(uint j=0;j < 37;++j)
			{
			for(uint k=0;k < ncolors;++k)
				{
					input_map(0,i,j,k) = 0.0;
				}
				
			}
			
		}
	
	for(uint i=0;i < indexedimage.size();++i)
		{
		for(uint j=0;j < indexedimage[i].second.size();++j)
			{
				input_map(0,indexedimage[i].first[0], indexedimage[i].first[1], j) = indexedimage[i].second[j];
			}	
		}

	//Actually run tensorflow
  	tensorflow::Status status = tfsession_->Run({ { "input_1", input_image }, {"input_2", input_b} },{ "k2tfout_0" }, {}, &outputs1);

	if (!status.ok()) 
		{
		std::cout << "Tensorflow Failed:" << std::endl;
  		std::cout << status.ToString() << "\n";
  		return ;
		}	
        float result_top = outputs1[0].flat<float>()(0);
        float result_qcd = outputs1[0].flat<float>()(1);

	itopdisc[jindex]=result_top/(result_top+result_qcd);
	jindex+=1;
  }


  //Add to jet userfloats
  auto outputs = std::make_unique<pat::JetCollection>();
  jindex=0;
  for (const auto &jet : *jets){
    pat::Jet newJet(jet);
    newJet.addUserFloat("Image"+extex_+":top", itopdisc[jindex]);
    outputs->push_back(newJet);
    jindex+=1;
  }

  // put into the event
  iEvent.put(std::move(outputs));
}

//define this as a plug-in
DEFINE_FWK_MODULE(ImageProducer);
