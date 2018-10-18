#ifndef BOOSTEDEVENTSHAPETAGGER_H
#define BOOSTEDEVENTSHAPETAGGER_H

// from: https://github.com/cms-ttbarAC/BESTAnalysis/blob/62a2fb6d514b45f7ae7b596a5bb2d41cca1bdaad/BoostedEventShapeTagger/interface/BoostedEventShapeTagger.h

#include <memory>
#include <fstream>

// CMS
#include "DataFormats/PatCandidates/interface/Jet.h"

// ROOT
#include "TLorentzVector.h"

// lwtnn
#include "lwtnn/LightweightNeuralNetwork.hh"
#include "lwtnn/parse_json.hh"


// BoostedEventShapeTagger Class
class BoostedEventShapeTagger {

  public:

    BoostedEventShapeTagger(const std::string& configFile, const std::string& dnnFile);
    ~BoostedEventShapeTagger();

    const std::map<std::string,double>& execute( const pat::Jet& jet );

    void getJetValues( const pat::Jet& jet );

    void pboost( TVector3 pbeam, TVector3 plab, TLorentzVector &pboo );

    void FWMoments( std::vector<TLorentzVector> particles, double (&outputs)[5] );

    float LegP(float x, int order);

    unsigned int getParticleID();

    void setConfigurations(const std::vector<std::string>& configurations);

    void read_file( const std::string &file_name, std::vector<std::string> &values, const std::string &comment="#" );

    bool str2bool( const std::string value );

  protected:

    // lwtnn
    std::unique_ptr<lwt::LightweightNeuralNetwork> m_lwtnn;
    std::map<std::string,double> m_BESTvars;
    std::map<std::string,double> m_NNresults;

    std::map<std::string,std::string> m_configurations; // map of configurations

    // kinematics
    float m_jetSoftDropMassMin; // [GeV] Jet soft drop mass minimum
    float m_jetPtMin;           // [GeV] Jet pT minimum
    unsigned int m_numSubjetsMin;    // minimum number of subjets
    unsigned int m_numDaughtersMin;  // minimum number of daughters

    // boosting to rest frames
    float m_radiusSmall;        // re-clustering jets
    float m_radiusLarge;        // re-clustering jets
    float m_reclusterJetPtMin;  // [GeV] re-clustering jet pT minimum

    float m_jetChargeKappa;     // weight for jet charge pT
    size_t m_maxJetSize;        // number of jets in re-clustering

    static constexpr float m_Wmass = 80.4;       // W mass [GeV]
    static constexpr float m_Zmass = 91.2;       // Z mass
    static constexpr float m_Hmass = 125.;       // Higgs mass
    static constexpr float m_Tmass = 172.5;      // Top mass

    std::map<std::string,std::string> m_defaultConfigs = {
             {"dnnFile",             "BESTAnalysis/BoostedEventShapeTagger/data/BEST_6bin_PUPPI.json"},
             {"radiusSmall",         "0.4"},
             {"radiusLarge",         "0.8"},
             {"reclusterJetPtMin",   "30.0"},
             {"jetSoftDropMassMin",  "40.0"},
             {"jetPtMin",            "500.0"},
             {"jetChargeKappa",      "0.6"},
             {"maxJetSize",          "4"},
             {"numSubjetsMin",       "2"},
             {"numDaughtersMin",     "2"} };
};

#endif
