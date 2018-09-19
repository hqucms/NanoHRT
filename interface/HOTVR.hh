//----------------------------------------------------------------------
//  
//  The Heavy Object Tagger with Variable R (HOTVR)
// 
//  This package combines variable R jet clustering with a 
//  mass jump veto. The resulting HOTVR jets have subjets, 
//  accessible through a helper class (HOTVRinfo).
//  Rejected clusters and jets without a mass jump 
//  can be accessed as well.
//
//  The code is based on the implementation of the ClusteringVetoPlugin 
//  version 1.0.0 (by Seng-Pei Liew and Martin Stoll)
//  and the VariableR plugin version 1.2.0 (by David Krohn, 
//  Gregory Soyez, Jesse Thaler and Lian-Tao Wang) in FastJet Contribs.
//  Please see the README file for more information.
//
//  For questions and comments, please contact: 
//     Tobias Lapsien  <tobias.lapsien@desy.de>
//     Roman Kogler    <roman.kogler@uni-hamburg.de>
//     Johannes Haller <johannes.haller@uni-hamburg.de>
//
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#ifndef __FASTJET_CONTRIB_HOTVR_HH__
#define __FASTJET_CONTRIB_HOTVR_HH__

#include <fastjet/internal/base.hh>
#include "HOTVRinfo.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include <fastjet/LimitedWarning.hh>

#include <queue>

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {

  //------------------------------------------------------------------------
  /// \class HOTVR
  /// 
  class HOTVR : public JetDefinition::Plugin {

  public:

    /// The strategy to be used with the clustering
    enum Strategy{
      Best,      ///< currently N2Tiled or N2Plain for FJ>3.2.0, NNH for FastJet<3.2.0
      N2Tiled,   ///< the default (faster in most cases) [requires FastJet>=3.2.0]
      N2Plain,   ///< [requires FastJet>=3.2.0]
      NNH,       ///< slower but already available for FastJet<3.2.0
    };

    // Result of veto function, from ClusteringVetoPlugin 1.0.0 
    enum VetoResult {
      CLUSTER,
      VETO,
      NOVETO
    };

    /// The clustering type is chosen by the "p" parameter
    /// generalised-kt algorithm. The definitions below are shorthand
    /// for the anti-kt, C/A and kt algorithm which also allow for
    /// backwards compatibility.
    static const double CALIKE;
    static const double KTLIKE;
    static const double AKTLIKE;

    /// Constructor that sets HOTVR algorithm parameters
    ///  - mu            mass jump threshold
    /// -  theta         mass jump strength
    ///  - min_r         minimum jet radius
    ///  - max_r         maximum jet radius
    ///  - rho           mass scale for effective radius (i.e. R ~ rho/pT)
    ///  - pt_sub        minimum pt of subjets
    ///  - clust_type    whether to use CA-like, kT-like, or anti-kT-like distance measure
    ///                  note that subjet finding has only been tested with CA-like clustering
    ///  - strategy      one of Best (the default), N2Tiled , N2Plain or NNH
    ///                  for FastJet>=3.2.0, the N2Tiled option is the default strategy, 
    ///                  for earlier FastJet versions NNH is used
    ///
    /// Example usage:
    /// \code
    /// HOTVR hotvr_plugin(mu, theta, min_r, max_r, rho, pt_sub, HOTVR::CALIKE);
    /// fastjet::JetDefinition jet_def(&hotvr_plugin);
    /// fastjet::ClusterSequence clust_seq(event, jet_def);
    /// vector<fastjet::PseudoJet> inclusive_jets = clust_seq.inclusive_jets(ptmin);
    /// std::vector<fastjet::PseudoJet> hotvr_jets;
    /// hotvr_jets = hotvr_plugin.get_jets();
    /// \endcode  
    HOTVR(double mu, double theta,double min_r, double max_r,double rho, double pt_sub, double clust_type, 
             Strategy requested_strategy = Best);
      
    // Virtual function from JetDefinition::Plugin that implements the algorithm
    void run_clustering(fastjet::ClusterSequence & cs) const;
   
    // Information string
    virtual string description() const;
      
    // NOTE: Required by JetDefinition::Plugin
    double R() const { return sqrt(_max_r2); }
    std::vector<PseudoJet>  get_jets(){return sorted_by_pt(_hotvr_jets);} //return the candidate jets
    std::vector<fastjet::PseudoJet> get_soft_cluster(){return _soft_cluster;} //return the clusters rejected by the comparison to a harder cluster
    std::vector<fastjet::PseudoJet> get_rejected_cluster(){return _rejected_cluster;} //return jets without massjumps
    std::vector<fastjet::PseudoJet> get_rejected_subjets(){return _rejected_subjets;} //return subjets rejected by the pT criterion
    void Reset(){ _hotvr_jets.clear(); _soft_cluster.clear();  _rejected_cluster.clear(); _subjets.clear(); _jets.clear(); _rejected_subjets.clear();}
      
  private:

    // bool for banner printout
    static bool _already_printed; 

    // Parameters of the HOTVR 
    double _mu, _theta, _min_r2, _max_r2, _max_r, _rho2, _pt_sub;
    double _clust_type;
    Strategy _requested_strategy;    
  
    // the jets and rejected clusters
    mutable std::vector<fastjet::PseudoJet> _jets;
    mutable std::vector<fastjet::PseudoJet> _hotvr_jets;
    mutable std::vector<fastjet::PseudoJet> _soft_cluster;
    mutable std::vector<fastjet::PseudoJet> _rejected_cluster;
    mutable std::vector<fastjet::PseudoJet> _rejected_subjets;
    mutable std::vector<std::vector<fastjet::PseudoJet> >  _subjets;

    // helper function to decide what strategy is best
    // the optimal strategy will depend on the multiplicity and _max_r
    Strategy best_strategy(unsigned int N) const;

    // implementation of the clustering using FastJet NN*** classes
    template<typename NN>
    void NN_clustering(ClusterSequence &cs, NN &nn) const;

    // veto condition
    VetoResult CheckVeto(const PseudoJet& j1, const PseudoJet& j2) const;
    
    // print welcome message with some information
    void print_banner();

  };
   
  //----------------------------------------------------------------------
  // classes to help run the HOTVR algorithm using NN-type classes
  // the code below is based on the implementation of the Variable R clustering
  // and has been taken from VariableR/VariableRPlugin.cc, version 1.2.1

  // class carrying particle-independent information
  class HOTVRNNInfo {
  public:
    HOTVRNNInfo(double rho2_in, double min_r2_in, double max_r2_in,
                    double clust_type_in)
      : _rho2(rho2_in), _min_r2(min_r2_in), _max_r2(max_r2_in),
        _clust_type(clust_type_in) {}
    
    double rho2()  const  {return _rho2; }
    double min_r2() const {return _min_r2;}
    double max_r2() const {return _max_r2;}
    double momentum_scale_of_pt2(double pt2) const {
      return pow(pt2,_clust_type);
    }
    
  private:
    double _rho2;   ///< constant (squared) that controls the overall R magnitude 
    double _min_r2; ///< minimal allowed radius squared
    double _max_r2; ///< maximal allowed radius squared
    double _clust_type; ///< cluster type (power "p" in distance measure)
  };

  // class carrying the minimal info required for the clustering
  class HOTVRBriefJet {

  public:

    void init(const PseudoJet & jet, HOTVRNNInfo *info) {
      _rap = jet.rap();
      _phi = jet.phi();
      double pt2 = jet.pt2();

      // get the effective "radius" Reff
      _beam_R2 = info->rho2()/pt2;
      if      (_beam_R2 > info->max_r2()){ _beam_R2 = info->max_r2();}
      else if (_beam_R2 < info->min_r2()){ _beam_R2 = info->min_r2();}

      // get the appropriate momentum scale
      _mom_factor2 = info->momentum_scale_of_pt2(pt2);
    }

    double geometrical_distance(const HOTVRBriefJet * jet) const {
      double dphi = std::abs(_phi - jet->_phi);
      double deta = (_rap - jet->_rap);
      if (dphi > pi) {dphi = twopi - dphi;}
      return dphi*dphi + deta*deta;
    }

    double geometrical_beam_distance() const { return _beam_R2; }

    double momentum_factor() const{ return _mom_factor2; }

    /// make this BJ class compatible with the use of NNH
    double distance(const HOTVRBriefJet * other_bj_jet){
      double mom1 = momentum_factor();
      double mom2 = other_bj_jet->momentum_factor();
      return (mom1<mom2 ? mom1 : mom2) * geometrical_distance(other_bj_jet);
    }
    double beam_distance(){
      return momentum_factor() * geometrical_beam_distance();
    }

    // the following are required by N2Tiled
    inline double rap() const{ return _rap;}
    inline double phi() const{ return _phi;}

  private:
    double _rap, _phi, _mom_factor2, _beam_R2;
  };
   
} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_HOTVR_HH__
