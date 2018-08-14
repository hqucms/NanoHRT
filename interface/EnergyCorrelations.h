/**
 * \file EnergyCorrelations.h
 * \brief Optimized code to calculate energy correlation functions. Based on code from fj-contrib
 * \author S.Narayanan
 */
#include "fastjet/PseudoJet.hh"
#include <vector>
#include <utility>
#include <tuple>
#include <map>
#include <type_traits>

#include "TMath.h"
#include "TString.h"

#ifndef PANDA_ECF_H
#define PANDA_ECF_H



/**
 * \brief delta-r-squared metric between two pseudojets
 * @param  j1 first jet
 * @param  j2 second jet
 * @return    \f$dR^2\f$
 */
 double jetDeltaR2(const fastjet::PseudoJet& j1, const fastjet::PseudoJet& j2);


/**
 * \brief in-place insertion sort given forward iterators
 */
template <typename FI>
void insertion_sort(FI b, FI e)
{
  int n = std::distance(b, e);
  for (int i = 1; i != n; ++i) {
    typename std::iterator_traits<FI>::value_type pivot = *(b+i); // current element is the pivot, make a copy
    int j = i - 1;
    while (j != -1 && pivot < *(b+j)) { // starting at pivot and going backwards
      *(b+j+1) = *(b+j); // if v[j] is less than pivot, move it up one index
      --j;
    }
    *(b+j+1) = pivot; // end up here either because j = -1 (pivot is smallest element) or
                      // because we found a j such that v[j]<pivot, in which case v[j+1] should be pivot
  }
}

template <typename T>
void insertion_sort(std::vector<T>& v) { insertion_sort(v.begin(), v.end()); }

template <typename FI>
typename std::iterator_traits<FI>::value_type min_e(FI b, FI e) { 
  int n = std::distance(b, e);
  typename std::iterator_traits<FI>::value_type r = *b;
  for (int i = 1; i != n; ++i) {
    r = std::min(r, *(b+i));
  }
  return r;
}

/**
 * \brief N-smallest elements
 * Find the N smallest elements and put them at the beginning of the vector
 * in smallest-to-biggest order
 * This works better than partial sorts when nmax << v.size()
 */
template <typename T>
void n_smallest(const std::vector<T>& v, int nmax) {
  auto b = v.begin(); auto e = v.end();
  for (int i = 0; i != nmax; ++i) {
    auto min_iter = std::min_element(b + i, e);
    std::iter_swap(b+i, min_iter); 
  }
}

/**
 * \brief in-place partial insertion sort
 * NB: this algorithm mangles the original vector above nmax!
 * This is done to minimize the number of comparison operations
 */
template <typename FI>
void partial_insertion_sort(FI b, FI e, int nmax)
{
  // clone of above, just don't do every comparison.
  int n = std::distance(b, e);
  for (int i = 1; i != n; ++i) {
    typename std::iterator_traits<FI>::value_type pivot = *(b+i); // make a copy
    int j = std::min(i - 1, nmax - 1); // whatever happens above nmax is junk
    while (j != -1 && pivot < *(b+j)) { 
      *(b+j+1) = *(b+j); 
      --j;
    }
    *(b+j+1) = pivot; 
  }
}

template <typename T>
void partial_insertion_sort(std::vector<T>& v, int nmax) { 
  partial_insertion_sort(v.begin(), v.end(), nmax);
}


class ECFCalculator {
public:
  enum param {
    oP=0, nP, bP, ecfP
  };
  typedef std::tuple<int, int, int, double> data_type;
  typedef std::tuple<int, int, int> pos_type;

  ECFCalculator(int maxN = 4,
             std::vector<float> bs = {0.5, 1, 2, 4});
  ~ECFCalculator() { }

  data_type access(int pos) const { return access(_oneToThree(pos)); }
  data_type access(pos_type pos) const;
  void calculate(const std::vector<fastjet::PseudoJet>&);

  // just a forward iterator
  class iterator {
  private:
    const ECFCalculator *_c;
    int _pos;
    data_type _data;

    void _access() { _data = _c->access(_pos); }
  public:
    iterator(const ECFCalculator *c, int pos = 0): _c(c), _pos(pos) { _access(); }
    iterator(const iterator& rhs): _c(rhs._c), _pos(rhs._pos), _data(rhs._data) { }
    ~iterator() { }

    iterator& operator++() { _pos++; _access(); return *this; }
    iterator operator++(int) { auto old(*this); ++(*this); return old; }
    iterator operator+(int n) const { return iterator(_c, _pos+n); }
    iterator& operator+=(int n) { _pos += n; _access(); return *this; }
    int operator-(const iterator& rhs) const { return this->_pos - rhs._pos; }
    bool operator==(const iterator& rhs) const { return this->_pos == rhs._pos; }
    bool operator!=(const iterator& rhs) const { return !( (*this) == rhs ); }
    const data_type& operator->() const { return _data; }
    template <int I>
      const typename std::tuple_element<I, data_type>::type& get() const { return std::get<I>(_data); }
  };

  iterator begin() const { return iterator(this, 0); }
  iterator end() const { return iterator(this, _bN * _nN * _oN); } 

private:
  void _set(pos_type pos, double x) { _ecfs[_threeToOne(pos)] = x; }
  pos_type _oneToThree(int pos) const;
  int _threeToOne(pos_type pos) const { return std::get<oP>(pos) 
                                               + (_oN * std::get<nP>(pos)) 
                                               + (_oN * _nN * std::get<bP>(pos)); }

  std::vector<float> _bs;
  std::vector<int> _ns, _os;
  const int _bN, _nN, _oN;
  std::vector<double> _ecfs;
  std::vector<double> pT; // these are member variables just to avoid
  std::vector<std::vector<double>> dR, dRBeta; // re-allocating memory 
};

#endif
