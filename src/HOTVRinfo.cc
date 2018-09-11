//----------------------------------------------------------------------
//
//  The Heavy Object Tagger with Variable R (HOTVR)
//
//  Helper class to access subjets and calculate jet observables.
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

#include "../interface/HOTVRinfo.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{
  double HOTVRinfo::mmin() const{ //calculate the minimum pairwise mass with the three leading subjets
  double mmin=0;
  double m12=0;
  double m01=0;
  double m02=0;
  if(_subjets.size()>2){
    m01 = 0;
    m01=(_subjets[0]+_subjets[1]).m();
    m02= 0;
    m02=(_subjets[0]+_subjets[2]).m();
    m12 = 0;
    m12 = (_subjets[1]+_subjets[2]).m();
  }
  mmin = std::min(m01,std::min(m02,m12));
  return mmin;
}


  double HOTVRinfo::ptfraction(int i) const { //calculate the pTfraction pT_{subjet,i}/pT_{jet}
  double ptfraction=0;
  if(_subjets.size()>i)
    ptfraction=_subjets.at(i).perp()/_parent.perp();
   return ptfraction;
}

  double HOTVRinfo::max_distance() const { //calculates the size of the jet, by finding the constituent with the largest distance. Can be used as an estimate for a jet radius
    std::vector<fastjet::PseudoJet> pfcands = _parent.constituents();
    double oldR=0;
    for(int j=0;j<pfcands.size();j++){
      double R;
      R=_parent.delta_R(pfcands[j]);
	    if(R>oldR){
	      oldR=R;
	    }
    }
    return oldR;
  }

}



FASTJET_END_NAMESPACE
