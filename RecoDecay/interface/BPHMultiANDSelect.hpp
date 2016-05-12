/*
 *  See header file for a description of this class.
 *
 *  $Date: 2015-07-06 18:40:19 $
 *  $Revision: 1.1 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//-----------------------
// This Class' Header --
//-----------------------


//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BPHAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "BPHAnalysis/RecoDecay/interface/BPHRecoSelect.h"
#include "BPHAnalysis/RecoDecay/interface/BPHMomentumSelect.h"
#include "BPHAnalysis/RecoDecay/interface/BPHVertexSelect.h"

//---------------
// C++ Headers --
//---------------


//-------------------
// Initializations --
//-------------------


//----------------
// Constructors --
//----------------
template<class T>
BPHMultiANDSelect<T>::BPHMultiANDSelect() {
}

//--------------
// Destructor --
//--------------
template<class T>
BPHMultiANDSelect<T>::~BPHMultiANDSelect() {
}

//--------------
// Operations --
//--------------
template<class T>
void BPHMultiANDSelect<T>::include( T& s, bool m ) {
  SelectElement e;
  e.selector = &s;
  e.mode     = m;
  selectList.push_back( e );
  return;
}


template<class T>
bool BPHMultiANDSelect<T>::accept( const reco::Candidate& cand,
                                   const BPHRecoBuilder* build ) const {
  return false;
}


template<class T>
bool BPHMultiANDSelect<T>::accept( const reco::Candidate& cand ) const {
  return false;
}


template<class T>
bool BPHMultiANDSelect<T>::accept( const BPHDecayMomentum& cand ) const {
  return false;
}


template<class T>
bool BPHMultiANDSelect<T>::accept( const BPHDecayVertex& cand ) const {
  return false;
}


template<>
bool BPHMultiANDSelect<BPHRecoSelect>::accept(
                                       const reco::Candidate& cand,
                                       const BPHRecoBuilder* build ) const;
template<>
bool BPHMultiANDSelect<BPHRecoSelect>::accept(
                                       const reco::Candidate& cand ) const;
template<>
bool BPHMultiANDSelect<BPHMomentumSelect>::accept(
                                           const BPHDecayMomentum& cand ) const;
template<>
bool BPHMultiANDSelect<BPHVertexSelect>::accept(
                                         const BPHDecayVertex& cand ) const;

