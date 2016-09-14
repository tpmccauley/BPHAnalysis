/*
 *  See header file for a description of this class.
 *
 *  \author Paolo Ronchese INFN Padova
 *
 */

//-----------------------
// This Class' Header --
//-----------------------
#include "BPHAnalysis/RecoDecay/interface/BPHMultiSelect.h"

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
// see interface/BPHMultiSelect.hpp

//--------------
// Destructor --
//--------------
// see interface/BPHMultiSelect.hpp

//--------------
// Operations --
//--------------
template<>
bool BPHMultiSelect<BPHRecoSelect>::accept(
                                    const reco::Candidate& cand,
                                    const BPHRecoBuilder* build ) const {
  return select( cand, build );
}


template<>
bool BPHMultiSelect<BPHRecoSelect>::accept(
                                    const reco::Candidate& cand ) const {
  return select( cand );
}


template<>
bool BPHMultiSelect<BPHMomentumSelect>::accept(
                                        const BPHDecayMomentum& cand ) const {
  return select( cand );
}


template<>
bool BPHMultiSelect<BPHVertexSelect>::accept(
                                      const BPHDecayVertex& cand ) const {
  return select( cand );
}

