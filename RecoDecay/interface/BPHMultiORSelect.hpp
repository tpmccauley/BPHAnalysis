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
//#include "BPHAnalysis/RecoDecay/interface/BPHMultiORSelect.h"

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
BPHMultiORSelect<T>::BPHMultiORSelect() {
}

//--------------
// Destructor --
//--------------
template<class T>
BPHMultiORSelect<T>::~BPHMultiORSelect() {
}

//--------------
// Operations --
//--------------
template<class T>
void BPHMultiORSelect<T>::include( T& s, bool m ) {
  SelectElement e;
  e.selector = &s;
  e.mode     = m;
  selectList.push_back( e );
  return;
}


template<class T>
bool BPHMultiORSelect<T>::accept( const reco::Candidate& cand ) const {
  return false;
}


template<class T>
bool BPHMultiORSelect<T>::accept( const reco::Candidate& cand,
                                  const BPHRecoBuilder* build ) const {
  return false;
}


template<class T>
bool BPHMultiORSelect<T>::accept( const BPHDecayMomentum& cand ) const {
  return false;
}


template<class T>
bool BPHMultiORSelect<T>::accept( const BPHDecayVertex& cand ) const {
  return false;
}


template<>
bool BPHMultiORSelect<BPHRecoSelect>::accept(
                                      const reco::Candidate& cand,
                                      const BPHRecoBuilder* build ) const;
template<>
bool BPHMultiORSelect<BPHRecoSelect>::accept(
                                      const reco::Candidate& cand ) const;
template<>
bool BPHMultiORSelect<BPHMomentumSelect>::accept(
                                          const BPHDecayMomentum& cand ) const;
template<>
bool BPHMultiORSelect<BPHVertexSelect>::accept(
                                        const BPHDecayVertex& cand ) const;

