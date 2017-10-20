/*
 *  See header file for a description of this class.
 *
 *  \author Paolo Ronchese INFN Padova
 *
 */

//-----------------------
// This Class' Header --
//-----------------------
#include "BPHAnalysis/SpecificDecay/interface/BPHBdToKxMuMuBuilder.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BPHAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "BPHAnalysis/RecoDecay/interface/BPHPlusMinusCandidate.h"
#include "BPHAnalysis/RecoDecay/interface/BPHRecoCandidate.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHMassSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHChi2Select.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHMassFitSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticleMasses.h"

//---------------
// C++ Headers --
//---------------
using namespace std;

//-------------------
// Initializations --
//-------------------


//----------------
// Constructors --
//----------------
BPHBdToKxMuMuBuilder::BPHBdToKxMuMuBuilder( const edm::EventSetup& es,
    const std::vector<BPHPlusMinusConstCandPtr>& oniaCollection,
    const std::vector<BPHPlusMinusConstCandPtr>&  kx0Collection ):
  oniaName( "Onia" ),
   kx0Name(  "Kx0" ),
  evSetup( &es ),
  jCollection( &oniaCollection ),
  kCollection( & kx0Collection ) {
  oniaSel = new BPHMassSelect   ( 1.00, 12.00 );
  mkx0Sel = new BPHMassSelect   ( 0.80, 1.00 );
  massSel = new BPHMassSelect   ( 3.50, 8.00 );
  chi2Sel = new BPHChi2Select   ( 0.02 );
  mFitSel = new BPHMassFitSelect( 4.00, 7.00 );
  massConstr = true;
  minPDiff = 1.0e-4;
  updated = false;
}

//--------------
// Destructor --
//--------------
BPHBdToKxMuMuBuilder::~BPHBdToKxMuMuBuilder() {
  delete oniaSel;
  delete mkx0Sel;
  delete massSel;
  delete chi2Sel;
  delete mFitSel;
}

//--------------
// Operations --
//--------------
vector<BPHRecoConstCandPtr> BPHBdToKxMuMuBuilder::build() {

  if ( updated ) return bdList;

  BPHRecoBuilder bBd( *evSetup );
  bBd.setMinPDiffererence( minPDiff );
  bBd.add( oniaName, *jCollection );
  bBd.add(  kx0Name, *kCollection );
  bBd.filter( oniaName, *oniaSel );
  bBd.filter(  kx0Name, *mkx0Sel );

  bBd.filter( *massSel );
  bBd.filter( *chi2Sel );
  if ( massConstr ) bBd.filter( *mFitSel );

  bdList = BPHRecoCandidate::build( bBd );
//
//  Apply kinematic constraint on the onia mass.
//  The operation is already performed when apply the mass selection,
//  so it's not repeated. The following code is left as example
//  for similar operations
//
//  int iBd;
//  int nBd = ( massConstr ? bdList.size() : 0 );
//  for ( iBd = 0; iBd < nBd; ++iBd ) {
//    BPHRecoCandidate* cptr( const_cast<BPHRecoCandidate*>(
//                            bdList[iBd].get() ) );
//    BPHRecoConstCandPtr onia = cptr->getComp( oniaName );
//    double oMass = onia->constrMass();
//    if ( oMass < 0 ) continue;
//    double sigma = onia->constrSigma();
//    cptr->kinematicTree( oniaName, oMass, sigma );
//  }
  updated = true;

  return bdList;

}

/// set cuts
void BPHBdToKxMuMuBuilder::setOniaMassMin( double m ) {
  updated = false;
  oniaSel->setMassMin( m );
  return;
}


void BPHBdToKxMuMuBuilder::setOniaMassMax( double m ) {
  updated = false;
  oniaSel->setMassMax( m );
  return;
}


void BPHBdToKxMuMuBuilder::setKxMassMin( double m ) {
  updated = false;
  mkx0Sel->setMassMin( m );
  return;
}


void BPHBdToKxMuMuBuilder::setKxMassMax( double m ) {
  updated = false;
  mkx0Sel->setMassMax( m );
  return;
}


void BPHBdToKxMuMuBuilder::setMassMin( double m ) {
  updated = false;
  massSel->setMassMin( m );
  return;
}


void BPHBdToKxMuMuBuilder::setMassMax( double m ) {
  updated = false;
  massSel->setMassMax( m );
  return;
}


void BPHBdToKxMuMuBuilder::setProbMin( double p ) {
  updated = false;
  chi2Sel->setProbMin( p );
  return;
}


void BPHBdToKxMuMuBuilder::setMassFitMin( double m ) {
  updated = false;
  mFitSel->setMassMin( m );
  return;
}


void BPHBdToKxMuMuBuilder::setMassFitMax( double m ) {
  updated = false;
  mFitSel->setMassMax( m );
  return;
}


void BPHBdToKxMuMuBuilder::setConstr( bool flag ) {
  updated = false;
  massConstr = flag;
  return;
}

/// get current cuts
double BPHBdToKxMuMuBuilder::getOniaMassMin() const {
  return oniaSel->getMassMin();
}


double BPHBdToKxMuMuBuilder::getOniaMassMax() const {
  return oniaSel->getMassMax();
}


double BPHBdToKxMuMuBuilder::getKxMassMin() const {
  return mkx0Sel->getMassMin();
}


double BPHBdToKxMuMuBuilder::getKxMassMax() const {
  return mkx0Sel->getMassMax();
}


double BPHBdToKxMuMuBuilder::getMassMin() const {
  return massSel->getMassMin();
}


double BPHBdToKxMuMuBuilder::getMassMax() const {
  return massSel->getMassMax();
}


double BPHBdToKxMuMuBuilder::getProbMin() const {
  return chi2Sel->getProbMin();
}


double BPHBdToKxMuMuBuilder::getMassFitMin() const {
  return mFitSel->getMassMin();
}


double BPHBdToKxMuMuBuilder::getMassFitMax() const {
  return mFitSel->getMassMax();
}


bool BPHBdToKxMuMuBuilder::getConstr() const {
  return massConstr;
}

