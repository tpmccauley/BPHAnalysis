/*
 *  See header file for a description of this class.
 *
 *  \author Paolo Ronchese INFN Padova
 *
 */

//-----------------------
// This Class' Header --
//-----------------------
#include "BPHAnalysis/SpecificDecay/interface/BPHBdToJPsiKsBuilder.h"

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
BPHBdToJPsiKsBuilder::BPHBdToJPsiKsBuilder( const edm::EventSetup& es,
    const std::vector<BPHPlusMinusConstCandPtr>& jpsiCollection,
    const std::vector<BPHPlusMinusConstCandPtr>&  k0sCollection ):
  jPsiName( "JPsi" ),
   k0sName(  "K0s" ),
  evSetup( &es ),
  jCollection( &jpsiCollection ),
  kCollection( & k0sCollection ) {
  jpsiSel = new BPHMassSelect   ( 2.80, 3.40 );
  mk0sSel = new BPHMassFitSelect( 0.00, 2.00 );
  massSel = new BPHMassSelect   ( 3.50, 8.00 );
  probMin = 0.02;
  mFitSel = new BPHMassFitSelect( jPsiName,
                                  BPHParticleMasses::jPsiMass,
                                  BPHParticleMasses::jPsiMWidth,
                                  5.00, 6.00 );	
  massConstr = true;
  minPDiff = 1.0e-4;
  updated = false;
}

//--------------
// Destructor --
//--------------
BPHBdToJPsiKsBuilder::~BPHBdToJPsiKsBuilder() {
  delete jpsiSel;
  delete mk0sSel;
  delete massSel;
  delete mFitSel;
}

//--------------
// Operations --
//--------------
vector<BPHRecoConstCandPtr> BPHBdToJPsiKsBuilder::build() {

  if ( updated ) return bdList;

  bdList.clear();
  dMap.clear();

  BPHRecoBuilder bBd( *evSetup );
  bBd.setMinPDiffererence( minPDiff );
  bBd.add( jPsiName, *jCollection );
  bBd.add(  k0sName, *kCollection );

  bBd.filter( jPsiName, *jpsiSel );
  bBd.filter(  k0sName, *mk0sSel );

  bBd.filter( *massSel );

  vector<BPHRecoConstCandPtr> tmpList = BPHRecoCandidate::build( bBd );
//
//  Apply kinematic constraint on the JPsi mass.
//
  int iBd;
  int nBd = tmpList.size();
  bdList.reserve( nBd );
  for ( iBd = 0; iBd < nBd; ++iBd ) {
    BPHRecoConstCandPtr ctmp = tmpList[iBd];
    BPHRecoCandidate* cptr = const_cast<BPHRecoCandidate*>( ctmp->clone() );
    BPHRecoConstCandPtr cand( cptr );
    cptr->setIndependentFit( k0sName );    // fit for K0s reconstruction
                                           // indipendent from other particles
    BPHRecoCandidate* kptr = const_cast<BPHRecoCandidate*>(
                      cptr->getComp( k0sName ).get() );
    kptr->setConstraint( BPHParticleMasses::k0sMass,
                         BPHParticleMasses::k0sMSigma );
    cptr->resetKinematicFit();
    if ( !mFitSel->accept( *cptr ) ) continue;
    const RefCountedKinematicVertex tdv = cptr->topDecayVertex();
    if ( tdv.get() == 0 ) continue;
    if ( !tdv->vertexIsValid() ) continue;
    reco::Vertex vtx( *tdv );
    if ( TMath::Prob( vtx.chi2(), lround( vtx.ndof() ) ) < probMin ) continue;
    dMap[cand->getComp( jPsiName ).get()] = ctmp->getComp( jPsiName ).get();
    dMap[cand->getComp(  k0sName ).get()] = ctmp->getComp(  k0sName ).get();
    bdList.push_back( cand );
  }
  updated = true;

  return bdList;

}

/// get original daughters map
const std::map<const BPHRecoCandidate*,
               const BPHRecoCandidate*>&
               BPHBdToJPsiKsBuilder::daughMap() const {
  return dMap;
}

/// set cuts
void BPHBdToJPsiKsBuilder::setJPsiMassMin( double m ) {
  updated = false;
  jpsiSel->setMassMin( m );
  return;
}


void BPHBdToJPsiKsBuilder::setJPsiMassMax( double m ) {
  updated = false;
  jpsiSel->setMassMax( m );
  return;
}


void BPHBdToJPsiKsBuilder::setK0MassMin( double m ) {
  updated = false;
  mk0sSel->setMassMin( m );
  return;
}


void BPHBdToJPsiKsBuilder::setK0MassMax( double m ) {
  updated = false;
  mk0sSel->setMassMax( m );
  return;
}


void BPHBdToJPsiKsBuilder::setMassMin( double m ) {
  updated = false;
  massSel->setMassMin( m );
  return;
}


void BPHBdToJPsiKsBuilder::setMassMax( double m ) {
  updated = false;
  massSel->setMassMax( m );
  return;
}


void BPHBdToJPsiKsBuilder::setProbMin( double p ) {
  updated = false;
  probMin = p;
  return;
}


void BPHBdToJPsiKsBuilder::setMassFitMin( double m ) {
  updated = false;
  mFitSel->setMassMin( m );
  return;
}


void BPHBdToJPsiKsBuilder::setMassFitMax( double m ) {
  updated = false;
  mFitSel->setMassMax( m );
  return;
}


void BPHBdToJPsiKsBuilder::setConstr( bool flag ) {
  updated = false;
  massConstr = flag;
  return;
}

/// get current cuts
double BPHBdToJPsiKsBuilder::getJPsiMassMin() const {
  return jpsiSel->getMassMin();
}


double BPHBdToJPsiKsBuilder::getJPsiMassMax() const {
  return jpsiSel->getMassMax();
}


double BPHBdToJPsiKsBuilder::getK0MassMin() const {
  return mk0sSel->getMassMin();
}


double BPHBdToJPsiKsBuilder::getK0MassMax() const {
  return mk0sSel->getMassMax();
}


double BPHBdToJPsiKsBuilder::getMassMin() const {
  return massSel->getMassMin();
}


double BPHBdToJPsiKsBuilder::getMassMax() const {
  return massSel->getMassMax();
}


double BPHBdToJPsiKsBuilder::getProbMin() const {
  return probMin;
}


double BPHBdToJPsiKsBuilder::getMassFitMin() const {
  return mFitSel->getMassMin();
}


double BPHBdToJPsiKsBuilder::getMassFitMax() const {
  return mFitSel->getMassMax();
}


bool BPHBdToJPsiKsBuilder::getConstr() const {
  return massConstr;
}

