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
std::set<const BPHRecoCandidate*> BPHRecoCandidate::allCand;

//----------------
// Constructors --
//----------------
BPHRecoCandidate::BPHRecoCandidate( const edm::EventSetup* es ):
  BPHDecayVertex( es ) {
  addCand();
}


BPHRecoCandidate::BPHRecoCandidate( const edm::EventSetup* es,
                  const BPHRecoBuilder::ComponentSet& compList ):
  BPHDecayMomentum( compList.daugMap, compList.compMap ),
  BPHDecayVertex( es ),
  BPHKinematicFit( BPHDecayMomentum::componentList() ) {
  addCand();

}

//--------------
// Destructor --
//--------------
BPHRecoCandidate::~BPHRecoCandidate() {
  // remove this object form the bookkeeping list
  allCand.erase( this );
}


//--------------
// Operations --
//--------------
std::vector<const BPHRecoCandidate*> BPHRecoCandidate::build(
                                     const BPHRecoBuilder& builder,
                                     double mass, double msig ) {
  // create a list of pointers to BPHRecoCandidate and fill it
  // with particle combinations selected by the BPHRecoBuilder
  std::vector<const BPHRecoCandidate*> cList;
  fill( cList, builder, mass, msig );
  return cList;
}


void BPHRecoCandidate::clear() {
  // repeatedly delete all objects
  // each objects is automatiaclly removed from the list when deleted
  while ( allCand.size() ) delete *allCand.begin();
  return;
}


void BPHRecoCandidate::addCand() {
  allCand.insert( this );
  return;
}

