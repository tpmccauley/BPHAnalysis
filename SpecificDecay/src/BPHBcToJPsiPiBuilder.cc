/*
 *  See header file for a description of this class.
 *
 *  \author Paolo Ronchese INFN Padova
 *
 */

//-----------------------
// This Class' Header --
//-----------------------
#include "BPHAnalysis/SpecificDecay/interface/BPHBcToJPsiPiBuilder.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BPHAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "BPHAnalysis/RecoDecay/interface/BPHPlusMinusCandidate.h"
#include "BPHAnalysis/RecoDecay/interface/BPHRecoCandidate.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticleNeutralVeto.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticlePtSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticleEtaSelect.h"
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
BPHBcToJPsiPiBuilder::BPHBcToJPsiPiBuilder( const edm::EventSetup& es,
    const std::vector<BPHPlusMinusConstCandPtr>& jpsiCollection,
    const BPHRecoBuilder::BPHGenericCollection*  pionCollection ):
  jPsiName( "JPsi" ),
  pionName( "Pion" ),
  evSetup( &es ),
  jCollection( &jpsiCollection ),
  pCollection(  pionCollection )  {
  jpsiSel = new BPHMassSelect       ( 2.80, 3.40 );
   pnVeto = new BPHParticleNeutralVeto;
    ptSel = new BPHParticlePtSelect (  3.0 );
   etaSel = new BPHParticleEtaSelect( 10.0 );
  massSel = new BPHMassSelect       ( 4.00, 9.00 );
  chi2Sel = new BPHChi2Select       ( 0.02 );
  mFitSel = new BPHMassFitSelect    ( jPsiName,
                                      BPHParticleMasses::jPsiMass,
                                      BPHParticleMasses::jPsiMWidth,
                                      6.00, 7.00 );
  massConstr = true;
  minPDiff = 1.0e-4;
  updated = false;
}

//--------------
// Destructor --
//--------------
BPHBcToJPsiPiBuilder::~BPHBcToJPsiPiBuilder() {
  delete jpsiSel;
  delete  pnVeto;
  delete   ptSel;
  delete  etaSel;
  delete massSel;
  delete chi2Sel;
  delete mFitSel;
}

//--------------
// Operations --
//--------------
vector<BPHRecoConstCandPtr> BPHBcToJPsiPiBuilder::build() {

  if ( updated ) return bcList;

  bcList.clear();

  BPHRecoBuilder bBc( *evSetup );
  bBc.setMinPDiffererence( minPDiff );
  bBc.add( jPsiName, *jCollection );
  bBc.add( pionName,  pCollection, BPHParticleMasses::pionMass,
                                   BPHParticleMasses::pionMSigma );
  bBc.filter( jPsiName, *jpsiSel );
  bBc.filter( pionName, * pnVeto );
  bBc.filter( pionName, *  ptSel );
  bBc.filter( pionName, * etaSel );

  bBc.filter( *massSel );
  if ( chi2Sel != 0 )
  bBc.filter( *chi2Sel );
  if ( massConstr ) bBc.filter( *mFitSel );

  bcList = BPHRecoCandidate::build( bBc );
//
//  Apply kinematic constraint on the JPsi mass.
//  The operation is already performed when apply the mass selection,
//  so it's not repeated. The following code is left as example
//  for similar operations
//
//  int iBc;
//  int nBc = ( massConstr ? bcList.size() : 0 );
//  for ( iBc = 0; iBc < nBc; ++iBc ) {
//    BPHRecoCandidate* cptr = const_cast<BPHRecoCandidate*>(
//                             bcList[iBc].get() );
//    BPHRecoConstCandPtr jpsi = cptr->getComp( jPsiName );
//    double jMass = jpsi->constrMass();
//    if ( jMass < 0 ) continue;
//    double sigma = jpsi->constrSigma();
//    cptr->kinematicTree( jPsiName, jMass, sigma );
//  }
  updated = true;
  return bcList;

}

/// set cuts
void BPHBcToJPsiPiBuilder::setJPsiMassMin( double m ) {
  updated = false;
  jpsiSel->setMassMin( m );
  return;
}


void BPHBcToJPsiPiBuilder::setJPsiMassMax( double m ) {
  updated = false;
  jpsiSel->setMassMax( m );
  return;
}


void BPHBcToJPsiPiBuilder::setPiPtMin( double pt ) {
  updated = false;
  ptSel->setPtMin( pt );
  return;
}


void BPHBcToJPsiPiBuilder::setPiEtaMax( double eta ) {
  updated = false;
  etaSel->setEtaMax( eta );
  return;
}


void BPHBcToJPsiPiBuilder::setMassMin( double m ) {
  updated = false;
  massSel->setMassMin( m );
  return;
}


void BPHBcToJPsiPiBuilder::setMassMax( double m ) {
  updated = false;
  massSel->setMassMax( m );
  return;
}


void BPHBcToJPsiPiBuilder::setProbMin( double p ) {
  updated = false;
  delete chi2Sel;
  chi2Sel = ( p < 0.0 ? 0 : new BPHChi2Select( p ) );
  return;
}


void BPHBcToJPsiPiBuilder::setMassFitMin( double m ) {
  updated = false;
  mFitSel->setMassMin( m );
  return;
}


void BPHBcToJPsiPiBuilder::setMassFitMax( double m ) {
  updated = false;
  mFitSel->setMassMax( m );
  return;
}


void BPHBcToJPsiPiBuilder::setConstr( bool flag ) {
  updated = false;
  massConstr = flag;
  return;
}

/// get current cuts
double BPHBcToJPsiPiBuilder::getJPsiMassMin() const {
  return jpsiSel->getMassMin();
}


double BPHBcToJPsiPiBuilder::getJPsiMassMax() const {
  return jpsiSel->getMassMax();
}


double BPHBcToJPsiPiBuilder::getPiPtMin() const {
  return ptSel->getPtMin();
}


double BPHBcToJPsiPiBuilder::getPiEtaMax() const {
  return etaSel->getEtaMax();
}


double BPHBcToJPsiPiBuilder::getMassMin() const {
  return massSel->getMassMin();
}


double BPHBcToJPsiPiBuilder::getMassMax() const {
  return massSel->getMassMax();
}


double BPHBcToJPsiPiBuilder::getProbMin() const {
  return ( chi2Sel == 0 ? -1.0 : chi2Sel->getProbMin() );
}


double BPHBcToJPsiPiBuilder::getMassFitMin() const {
  return mFitSel->getMassMin();
}


double BPHBcToJPsiPiBuilder::getMassFitMax() const {
  return mFitSel->getMassMax();
}


bool BPHBcToJPsiPiBuilder::getConstr() const {
  return massConstr;
}

