/*
 *  See header file for a description of this class.
 *
 *  $Date: 2015-07-03 10:08:22 $
 *  $Revision: 1.1 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//-----------------------
// This Class' Header --
//-----------------------
#include "BPHAnalysis/RecoDecay/interface/BPHRecoCandidate.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BPHAnalysis/RecoDecay/interface/BPHRecoBuilder.h"

//---------------
// C++ Headers --
//---------------
#include <iostream>

using namespace std;

//-------------------
// Initializations --
//-------------------


//----------------
// Constructors --
//----------------
BPHRecoCandidate::BPHRecoCandidate( const edm::EventSetup* es ):
  BPHDecayVertex( es ) {
}


BPHRecoCandidate::BPHRecoCandidate( const edm::EventSetup* es,
                  const BPHRecoBuilder::ComponentSet& compList ):
  BPHDecayMomentum( compList.daugMap, compList.compMap ),
  BPHDecayVertex( this, es ),
  BPHKinematicFit( this ) {
}

//--------------
// Destructor --
//--------------
BPHRecoCandidate::~BPHRecoCandidate() {
}


//--------------
// Operations --
//--------------
std::vector<BPHRecoConstCandPtr> BPHRecoCandidate::build(
                                 const BPHRecoBuilder& builder,
                                 double mass, double msig ) {
  // create a list of pointers to BPHRecoCandidate and fill it
  // with particle combinations selected by the BPHRecoBuilder
  std::vector<BPHRecoConstCandPtr> cList;
  fill( cList, builder, mass, msig );
  return cList;
}

