/**
 * \file EnergyCorrelations.cc
 * \brief Re-optimized code to calculate energy correlation functions. 
 * \author S.Narayanan
 */
#include "../interface/EnergyCorrelations.h"
#define PI 3.141592654

#include <iostream>

using namespace std;
typedef ECFCalculator C;

double jetDeltaR2(const fastjet::PseudoJet& j1, const fastjet::PseudoJet& j2) 
{
    double dEta{j1.eta()-j2.eta()}; 
    double dPhi{j1.phi()-j2.phi()};

    if (dPhi<-PI)
        dPhi = 2*PI+dPhi;
    else if (dPhi>PI)
        dPhi = -2*PI+dPhi;

    return dEta*dEta + dPhi*dPhi;
}

C::ECFCalculator(int maxN, vector<float> bs):
  _bs(bs),
  _os({1,2,3}),
  _bN(_bs.size()),
  _nN(maxN),
  _oN(_os.size())
{
  _ns = vector<int>(maxN);
  for (int i = 0; i != maxN; ++i)
    _ns[i] = i + 1;
  _ecfs.resize(_nN * _oN * _bN);
}

C::data_type C::access(C::pos_type pos) const
{
  if (_threeToOne(pos) == _oN * _nN * _bN) {
    // don't throw error - fail silently, to allow end() to be defined
    return make_tuple(-1,-1,-1,-1);
  }
  return make_tuple(get<oP>(pos),
                    get<nP>(pos),
                    get<bP>(pos),
                    _ecfs.at(_threeToOne(pos)));
}

C::pos_type C::_oneToThree(int pos) const
{
  int oI = pos % _oN;
  pos -= oI; pos /= _oN;
  int nI = pos % _nN;
  pos -= nI; pos /= _nN;
  int bI = pos;
  return make_tuple(oI, nI, bI);
}

void C::calculate(const vector<fastjet::PseudoJet>& particles)
{
  if (_nN == 0 || _bN == 0 || _oN == 0)
    return;

  int nParticles = particles.size();

  // cache kinematics
  if (nParticles > (int)pT.size()) { 
    pT.resize(nParticles, 0); dR.resize(nParticles); dRBeta.resize(nParticles);
    for (int iP=0; iP!=nParticles; ++iP) {
      dR[iP].resize(nParticles, 0);
      dRBeta[iP].resize(nParticles, 0);
    }
  }

  for (int iP=0; iP!=nParticles; ++iP) {
    const fastjet::PseudoJet& pi = particles[iP];
    pT[iP] = pi.perp();
    for (int jP=0; jP!=iP; ++jP) {
      if (iP == jP) {
        dR[iP][jP] = 0;
      } else { 
        const fastjet::PseudoJet& pj = particles[jP];
        dR[iP][jP] = jetDeltaR2(pi,pj);
      }
    }
  }

  // get the normalization factor
  double baseNorm{0};
  for (auto& pt : pT)
    baseNorm += pt;
  double norm2{pow(baseNorm, 2)};
  double norm3{pow(baseNorm, 3)};
  double norm4{pow(baseNorm, 4)};

  vector<vector<double>> vals(_nN);
  for (int nI = 0; nI != _nN; ++nI) {
    vals[nI].resize(_oN, 0);
  }

  for (int bI = 0; bI != _bN; ++bI) {
    
    // reweight angles
    double halfBeta{_bs[bI] / 2.};
    for (int iP=0; iP!=nParticles; ++iP) {
      for (int jP=0; jP!=iP; ++jP) {
        if (iP == jP)
          dRBeta[iP][jP] = 0;
        else
          dRBeta[iP][jP] = pow(dR[iP][jP], halfBeta);
      }
    }

    // now we compute the ECFNs
    // trivial case, n = 1
    for (int oI = 0; oI != _oN; ++oI) {
      _set(make_tuple(oI, 0, bI), 1);
    }

    if (_nN < 2)
      continue;

    // reset the accumulators
    for (auto& v : vals) 
      std::fill(v.begin(), v.end(), 0);

    // now we loop
    vector<double> angles3(3), angles4(6); // N(N-1)/2
    for (int iP = 0; iP != nParticles; ++iP) {
      for (int jP = 0; jP != iP; ++jP) {
        const double pt_ij = pT[iP] * pT[jP];
        const double angle_ij = dRBeta[iP][jP];

        vals[1][0] += pt_ij * angle_ij;

        if (_nN > 2) { 
          for (int kP = 0; kP != jP; ++kP) {
            const double angle_ik = dRBeta[iP][kP], angle_jk = dRBeta[jP][kP];
            const double pt_ijk = pt_ij * pT[kP];

            angles3[0] = angle_ij; 
            angles3[1] = angle_ik; 
            angles3[2] = angle_jk;

            insertion_sort(angles3);
            
            // unrolling this appears to be faster than a for loop
            double inc3 = pt_ijk * angles3[0];
            vals[2][0] += inc3;
            inc3 *= angles3[1]; vals[2][1] += inc3;
            inc3 *= angles3[2]; vals[2][2] += inc3;


            if (_nN > 3) {
              // Set these outside of the  lP loop as long as the sorting
              // algorithm we use is not in-place 
              // If we are using an in-place sort, this needs to be moved down
              angles4[0] = angle_ij;
              angles4[1] = angle_ik;
              angles4[2] = angle_jk;

              for (int lP = 0; lP != kP; ++lP) {
                const double pt_ijkl = pt_ijk * pT[lP];
                angles4[3] = dRBeta[iP][lP]; 
                angles4[4] = dRBeta[jP][lP];
                angles4[5] = dRBeta[kP][lP];

                // The best-performing algorithm to find the 2- or 
                // 3-smallest elements appears to be as below. 
                // I have a number of other options in the header that were tested, 
                // all demonstrably worse. 
                // The exact performance of this sorting is critical, since it's called 
                // O(4e8) times/event.
                // "This is the worst sorting algorithm, except for all the other
                // ones that have been tried" -Winsort Churchill
                double angle1=999, angle2=999;
                int idx_1=-1;
                for (int iA = 0; iA != 6; ++iA) {
                  if (angles4[iA] < angle1) {
                    angle1 = angles4[iA];
                    idx_1 = iA;
                  }
                }
                for (int iA = 0; iA != 6; ++iA) {
                  if (iA == idx_1)
                    continue;
                  if (angles4[iA] < angle2) {
                    angle2 = angles4[iA];
                  }
                }


                // again prefer to unroll
                double inc4 = pt_ijkl * angle1;
                vals[3][0] += inc4;
                inc4 *= angle2; vals[3][1] += inc4;
                
              } // l
            } // N > 3

          } // k
        } // N > 2
      } // j
    } // i

  
    // set the values
    double val2 = vals[1][0];
    val2 /= norm2;
    for (int oI = 0; oI != _oN; ++oI) {
      _set(make_tuple(oI, 1, bI), val2);
      _set(make_tuple(oI, 2, bI), vals[2][oI] / norm3);
      _set(make_tuple(oI, 3, bI), vals[3][oI] / norm4);
    }
  } // beta loop
}

