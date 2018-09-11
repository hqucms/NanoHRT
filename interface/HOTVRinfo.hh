//----------------------------------------------------------------------
//  
//  The Heavy Object Tagger with Variable R (HOTVR)
// 
//  Helper class to access subjets and calculate jet observables.
//
//  For questions and comments, please contact: 
//      Tobias Lapsien  <tobias.lapsien@desy.de>
//      Roman Kogler    <roman.kogler@uni-hamburg.de>
//      Johannes Haller <johannes.haller@uni-hamburg.de>
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

#ifndef __FASTJET_CONTRIB_HOTVRINFO_HH__
#define __FASTJET_CONTRIB_HOTVRINFO_HH__
#include "fastjet/PseudoJet.hh"
FASTJET_BEGIN_NAMESPACE  

namespace contrib{

  class HOTVRinfo : public fastjet::PseudoJet::UserInfoBase {
  public:

    // inline constructor
    // construct the helper class with the subjets and the the parent jet itself
    HOTVRinfo(fastjet::PseudoJet jet, std::vector<PseudoJet> subjets): _subjets(subjets), _parent(jet) {}; 

    // size of the jet, defined by the constituent with the largest distance DeltaR from 
    // the jet axis, can can be used as an estimate for a jet radius
    double max_distance() const;

    // vector of the associated subjets
    std::vector<fastjet::PseudoJet> subjets()const {return _subjets;};

    // number of associated subjets
    double nsubjets() const {return _subjets.size();};

    // minimum pairwise mass of the three leading subjets
    double mmin() const;

    // ptfraction of subjet i, given by p_T,i / p_T,fatjet
    double ptfraction(int i) const;
       
  private:
    std::vector<fastjet::PseudoJet> _subjets;
    fastjet::PseudoJet _parent;
  };
  

}
FASTJET_END_NAMESPACE
#endif // __FASTJET_CONTRIB_HOTVRINFO_HH__
