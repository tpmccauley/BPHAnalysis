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
#include "BPHAnalysis/RecoDecay/interface/BPHMultiANDSelect.h"

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
// see interface/BPHMultiANDSelect.hpp

//--------------
// Destructor --
//--------------
// see interface/BPHMultiANDSelect.hpp

//--------------
// Operations --
//--------------
template<>
bool BPHMultiANDSelect<BPHRecoSelect>::accept(
                                       const reco::Candidate& cand,
                                       const BPHRecoBuilder* build ) const {
  int i;
  int n = selectList.size();
  for ( i = 0; i < n; ++i ) {
    const SelectElement& e = selectList[i];
    if ( e.selector->accept( cand, build ) != e.mode ) return false;
  }
  return true;
}


template<>
bool BPHMultiANDSelect<BPHRecoSelect>::accept(
                                       const reco::Candidate& cand ) const {
  int i;
  int n = selectList.size();
  for ( i = 0; i < n; ++i ) {
    const SelectElement& e = selectList[i];
    if ( e.selector->accept( cand ) != e.mode ) return false;
  }
  return true;
}


template<>
bool BPHMultiANDSelect<BPHMomentumSelect>::accept(
                                           const BPHDecayMomentum& cand ) const {
  int i;
  int n = selectList.size();
  for ( i = 0; i < n; ++i ) {
    const SelectElement& e = selectList[i];
    if ( e.selector->accept( cand ) != e.mode ) return false;
  }
  return true;
}


template<>
bool BPHMultiANDSelect<BPHVertexSelect>::accept(
                                         const BPHDecayVertex& cand ) const {
  int i;
  int n = selectList.size();
  for ( i = 0; i < n; ++i ) {
    const SelectElement& e = selectList[i];
    if ( e.selector->accept( cand ) != e.mode ) return false;
  }
  return true;
}

