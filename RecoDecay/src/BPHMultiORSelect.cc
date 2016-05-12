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
#include "BPHAnalysis/RecoDecay/interface/BPHMultiORSelect.h"

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
// see interface/BPHMultiORSelect.hpp

//--------------
// Destructor --
//--------------
// see interface/BPHMultiORSelect.hpp

//--------------
// Operations --
//--------------
template<>
bool BPHMultiORSelect<BPHRecoSelect>::accept(
                                      const reco::Candidate& cand,
                                      const BPHRecoBuilder* build ) const {
  int i;
  int n = selectList.size();
  for ( i = 0; i < n; ++i ) {
    const SelectElement& e = selectList[i];
    if ( e.selector->accept( cand, build ) == e.mode ) return true;
  }
  return false;
}


template<>
bool BPHMultiORSelect<BPHRecoSelect>::accept(
                                      const reco::Candidate& cand ) const {
  int i;
  int n = selectList.size();
  for ( i = 0; i < n; ++i ) {
    const SelectElement& e = selectList[i];
    if ( e.selector->accept( cand ) == e.mode ) return true;
  }
  return false;
}


template<>
bool BPHMultiORSelect<BPHMomentumSelect>::accept(
                                          const BPHDecayMomentum& cand ) const {
  int i;
  int n = selectList.size();
  for ( i = 0; i < n; ++i ) {
    const SelectElement& e = selectList[i];
    if ( e.selector->accept( cand ) == e.mode ) return true;
  }
  return false;
}


template<>
bool BPHMultiORSelect<BPHVertexSelect>::accept(
                                        const BPHDecayVertex& cand ) const {
  int i;
  int n = selectList.size();
  for ( i = 0; i < n; ++i ) {
    const SelectElement& e = selectList[i];
    if ( e.selector->accept( cand ) == e.mode ) return true;
  }
  return false;
}

