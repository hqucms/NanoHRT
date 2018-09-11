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

#include "../interface/HOTVR.hh"

#include <cstdio>
#include "math.h"
#include <iomanip>
#include <cmath>
#include <sstream>
#include <memory>

#include <fastjet/NNH.hh>
#if FASTJET_VERSION_NUMBER >= 30200
#include <fastjet/NNFJN2Plain.hh>
#include <fastjet/NNFJN2Tiled.hh>
#endif


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {

  // make sure that welcome banner is printed only once
  bool HOTVR::_already_printed = false;

  /// The clustering type is chosen by the "p" parameter
  /// generalised-kt algorithm. The definitions below are shorthand
  /// for the anti-kt, C/A and kt algorithm which also allow for
  /// backwards compatibility.
  const double HOTVR::CALIKE  =  0.0;
  const double HOTVR::KTLIKE  =  1.0;
  const double HOTVR::AKTLIKE = -1.0;

  // Constructor that sets HOTVR algorithm parameters
  //  - mu            mass jump threshold
  // -  theta         mass jump strength
  //  - min_r         minimum jet radius
  //  - max_r         maximum jet radius
  //  - rho           mass scale for effective radius (i.e. R ~ rho/pT)
  //  - pt_sub        minimum pt of subjets
  //  - clust_type    whether to use CA-like, kT-like, or anti-kT-like distance measure
  //                  note that subjet finding has only been tested with CA-like clustering
  //  - strategy      one of Best (the default), N2Tiled , N2Plain or NNH
  HOTVR::HOTVR(double mu, double theta, double min_r, double max_r, double rho, double pt_sub, double clust_type, Strategy requested_strategy)
    :_mu(mu), _theta(theta), _min_r2(min_r*min_r), _max_r2(max_r*max_r), _max_r(max_r),
     _rho2(rho*rho), _pt_sub(pt_sub),
     _clust_type(clust_type), _requested_strategy(requested_strategy)
  {
    if (!_already_printed){
//      print_banner();
      _already_printed = true;
    }

    // some sanity checks
    if (mu < 0.0) throw Error("HOTVR: mu must be positive.");
    if (theta > 1.0 || theta < 0.0) throw Error("HOTVR: theta must be in [0.0,1.0].");
    if (max_r < 0.0) throw Error("HOTVR: Maximum radius must be positive.");
    if (min_r < 0.0) throw Error("HOTVR: Minimum radius must be positive.");
    if (min_r>max_r)  throw Error("HOTVR: Minimum radius must be smaller than maximum radius.");
    if (rho<0)  throw Error("HOTVR: Rho must be positive.");
    if (pt_sub<0)  throw Error("HOTVR: pT threshold must be positive.");

    // decide the strategy
#if FASTJET_VERSION_NUMBER < 30200
    // this is only supported for the Best and Native strategies
    if ((requested_strategy!=Best) && (requested_strategy!=NNH))
      throw Error("HOTVR: with FastJet<3.2.0, Best and NNH are the only supported strategies.");
#endif

  }

  void HOTVR::run_clustering(ClusterSequence & cs) const {

    // set up NNH
    HOTVRNNInfo nninfo(_rho2, _min_r2, _max_r2, _clust_type);

    // the following code has been written by G. Soyez and is taken from
    // VariableR/VariableRPlugin.cc, version 1.2.1
    // -> make use of the NN-type clustering in FastJet 3.2 and higher
#if FASTJET_VERSION_NUMBER >= 30200

    // set up clustering strategy
    Strategy strategy = _requested_strategy;

    // decide the best option upon request
    if (_requested_strategy==Best){
      strategy = best_strategy(cs.jets().size());
    }

    if (strategy==N2Tiled){
      NNFJN2Tiled<HOTVRBriefJet,HOTVRNNInfo> nnt(cs.jets(), _max_r, &nninfo);
      NN_clustering(cs, nnt);
    } else if (strategy==N2Plain){
      NNFJN2Plain<HOTVRBriefJet,HOTVRNNInfo> nnp(cs.jets(), &nninfo);
      NN_clustering(cs, nnp);
    } else { // NNH is the only option left
#endif
      fastjet::NNH<HOTVRBriefJet,HOTVRNNInfo> nnh(cs.jets(), &nninfo);
      NN_clustering(cs, nnh);
#if FASTJET_VERSION_NUMBER >= 30200
    }
#endif
  }



  //---------------------------------------------------------------------
  // the actual implementation of the clustering using the NN helpers
  template<typename NN>
  void HOTVR::NN_clustering(ClusterSequence &cs, NN &nn) const{

    // loop over pseudojets
    int njets = cs.jets().size();
    while (njets > 0) {
      bool existing_i=false;
      bool existing_j=false;
      int i(-1), j(-1);
      double dij = nn.dij_min(i, j);

      if(j < 0) {//diB is smallest
	      cs.plugin_record_iB_recombination(i,dij);
	      _jets.push_back(cs.jets()[i]);
	      _jets[_jets.size()-1].set_user_index(-2);

	      bool set=false;
	      std::vector<fastjet::PseudoJet> subjets;
	      for(unsigned o=0;o<_jets.size();o++){
	        if(_jets[o].user_index()==i) {
	          _jets[_jets.size()-1].set_user_index(i*100);
	          if(!set) {
	            _hotvr_jets.push_back(cs.jets()[i]);//save jet candidates (jets with massjumps)
	            _hotvr_jets[_hotvr_jets.size()-1].set_user_index(_hotvr_jets.size()-1);
	            set=true;
	          }
	          subjets.push_back(_jets[o]);//save subjets
	        }
        }

        if(_jets[_jets.size()-1].user_index()!=-2) {
	        _subjets.push_back(sorted_by_pt(subjets));
	        _hotvr_jets.at(_hotvr_jets.size()-1).set_user_info(
          new HOTVRinfo(_hotvr_jets.at(_hotvr_jets.size()-1),
                        _subjets[_hotvr_jets.at(_hotvr_jets.size()-1).user_index()]) );
	      } else {
          _rejected_cluster.push_back(cs.jets()[i]); //fill vector with jets that have no massjumps
        }

        nn.remove_jet(i);//remove jet from list
	      njets--;
	     continue;
      }


      int k=-1;//dij is smallest
      switch ( CheckVeto ( cs.jets()[i],cs.jets()[j] ) ) {//below the massjump threshold mu: cluster

      case CLUSTER: {
	      k=-1;
	      cs.plugin_record_ij_recombination(i, j, dij, k);//cluster i and j
	      nn.merge_jets(i, j, cs.jets()[k], k);
	      njets--;
	      break;
      }

      case VETO: { //above massjump threshold mu: massjump found
	      k=-1;
	  	  if(cs.jets()[i].pt()<_pt_sub) {
	        cs.plugin_record_iB_recombination(i,dij);     //subjet i below pT threshold?
		      _rejected_subjets.push_back(cs.jets().at(i)); //save rejected subjets here
	        nn.remove_jet(i);
	        njets--;
	      }
	      if(cs.jets()[j].pt()<_pt_sub) {
	        cs.plugin_record_iB_recombination(j,dij);     //subjet j below pT threshold?
		      _rejected_subjets.push_back(cs.jets().at(j)); //save rejected subjets here
	        nn.remove_jet(j);
	        njets--;
	      }
	      if(cs.jets()[i].pt()>=_pt_sub && cs.jets()[j].pt()>=_pt_sub) {//check if the subjet pT is higher than the threshold
	        cs.plugin_record_ij_recombination(i, j, dij, k);

  	      for(unsigned o=0;o<_jets.size();o++){
	          if(_jets[o].user_index()==j) {
	            _jets[o].set_user_index(k);
	            existing_j=true;
	          }
	          if( _jets[o].user_index()==i){
	            _jets[o].set_user_index(k);
	            existing_i=true;
	          }
	        }
	        if(!existing_j){
	          _jets.push_back(cs.jets()[j]);
	          _jets[_jets.size()-1].set_user_index(k);
	        }
	        if(!existing_i){
	          _jets.push_back(cs.jets()[i]);
	          _jets[_jets.size()-1].set_user_index(k);
	        }
	        nn.merge_jets(i, j, cs.jets()[k], k);
	        njets--;
	        }
	      break;
      }

      case NOVETO: {  //above massjump threshold mu:  no massjump found
	      if(cs.jets()[i].m()<cs.jets()[j].m()) {
	        cs.plugin_record_iB_recombination(i,dij);
	        _soft_cluster.push_back(cs.jets()[i]); //fill vector with jets that were rejected
	        nn.remove_jet(i);
	      } else {
	        cs.plugin_record_iB_recombination(j,dij);
	        _soft_cluster.push_back(cs.jets()[j]); //fill vector with jets that were rejected
	        nn.remove_jet(j);
	      }
	      njets--;
      }

      } // end switch over veto condition

    } // end while loop over pseudojets

  }


  string HOTVR::description() const{
    stringstream sstream("");

    sstream << "HOTVR (1606.04961), ";

    if (_clust_type<0) {
      sstream << "AKT";
    } else if (_clust_type==0.) {
       sstream << "CA";
    } else {
      sstream << "KT";
    }
    sstream << "-like";
    sstream << fixed << setprecision(1) << ", theta=" << _theta;
    sstream << ", mu=" << _mu;
    sstream << ", max_r=" << sqrt(_max_r2);
    sstream << ", min_r=" << sqrt(_min_r2);
    sstream << ", rho=" << sqrt(_rho2);
    sstream << ", pt_sub=" << _pt_sub;

    switch (_requested_strategy){
    case Best:    sstream << ", strategy=Best"; break;
    case N2Tiled: sstream << ", strategy=N2Tiled"; break;
    case N2Plain: sstream << ", strategy=N2Plain"; break;
    case NNH:     sstream << ", strategy=NNH"; break;
    };

    return sstream.str();
  }


  void HOTVR::print_banner(){
    std::cout << "#------------------------------------------------------------------------------\n";
    std::cout << "#                                 HOTVR 1.0.0                                  \n";
    std::cout << "#                      T. Lapsien, R. Kogler, J. Haller                        \n";
    std::cout << "#                                arXiv:1606.04961                              \n";
    std::cout << "#                                                                              \n";
    std::cout << "# The Heavy Object Tagger with Variable R.                                     \n";
    std::cout << "#                                                                              \n";
    std::cout << "# If you use this package, please cite the following papers:                   \n";
    std::cout << "#                                                                              \n";
    std::cout << "# - Tobias Lapsien, Roman Kogler, Johannes Haller,                             \n";
    std::cout << "#   \"A new tagger for hadronically decaying heavy particles at the LHC\",     \n";
    std::cout << "#   arXiv:1606.04961                                                           \n";
    std::cout << "#                                                                              \n";
    std::cout << "# - David Krohn, Jesse Thaler, Lian-Tao Wang, \"Jets with Variable R\",        \n";
    std::cout << "#   JHEP 06, 059 (2009) [arXiv:0903.0392 [hep-ph]]                             \n";
    std::cout << "#                                                                              \n";
    std::cout << "# - Martin Stoll, \"Vetoed jet clustering: The mass-jump algorithm\",          \n";
    std::cout << "#   JHEP 04, 111 (2015) [arXiv:1410.4637 [hep-ph]]                             \n";
    std::cout << "#                                                                              \n";
    std::cout << "# The HOTVR algorithm is provided without warranty under the terms of the      \n";
    std::cout << "# GNU GPLv2. See COPYING file for details.                                     \n";
    std::cout << "#------------------------------------------------------------------------------\n";
}

  // check the mass jump terminating veto
  // this function is taken from ClusteringVetoPlugin 1.0.0
  HOTVR::VetoResult HOTVR::CheckVeto(const PseudoJet& j1, const PseudoJet& j2) const {

    PseudoJet combj = j1+j2;

    double mj1 = abs (j1.m());
    double mj2 = abs (j2.m());
    double mcombj = abs (combj.m());

    if (mcombj < _mu)  return CLUSTER; // recombine
    else if (_theta*mcombj > max(mj1,mj2)) return VETO;
    else return NOVETO; // mass jump

  }

  // decide the optimal strategy
  // taken from VariableR/VariableRPlugin.cc, version 1.2.1
  HOTVR::Strategy HOTVR::best_strategy(unsigned int N) const{
#if FASTJET_VERSION_NUMBER >= 30200
    // use the FastJet (v>/3.1) transition between N2Plain and N2Tiled
    if (N <= 30 || N <= 39.0/(max(_max_r, 0.1) + 0.6)) return N2Plain;
    else return N2Tiled;
#else
    return NNH;
#endif
  }


} // namespace contrib

FASTJET_END_NAMESPACE
