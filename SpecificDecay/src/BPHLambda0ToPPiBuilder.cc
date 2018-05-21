/*
 *  See header file for a description of this class.
 *
 *  \author Paolo Ronchese INFN Padova
 *
 */

//-----------------------
// This Class' Header --
//-----------------------
#include "BPHAnalysis/SpecificDecay/interface/BPHLambda0ToPPiBuilder.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BPHAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "BPHAnalysis/RecoDecay/interface/BPHPlusMinusCandidate.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticlePtSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticleEtaSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHMassSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHMassSymSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHChi2Select.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticleMasses.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"

//---------------
// C++ Headers --
//---------------
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string>
using namespace std;

//-------------------
// Initializations --
//-------------------


//----------------
// Constructors --
//----------------

BPHLambda0ToPPiBuilder::BPHLambda0ToPPiBuilder(
               const edm::EventSetup& es,
               const BPHRecoBuilder::BPHGenericCollection* protonCollection,
               const BPHRecoBuilder::BPHGenericCollection*   pionCollection ):
  protonName( "Proton" ),
    pionName(   "Pion" ),
  evSetup( &es ),
  prCollection( protonCollection ),
  piCollection(   pionCollection ), 
  vCollection( 0 ),
  rCollection( 0 ),
  sList( "" ) {
    ptSel = new BPHParticlePtSelect (  0.7 );
   etaSel = new BPHParticleEtaSelect( 10.0 );
  massSel = new BPHMassSelect( 0.80, 1.40 ) ;
  chi2Sel = new BPHChi2Select( 0.0 );
  updated = false;
}

BPHLambda0ToPPiBuilder::BPHLambda0ToPPiBuilder(
               const edm::EventSetup& es,
               const std::vector<reco::VertexCompositeCandidate>*
                     v0Collection,
               const std::string& searchList ):
  protonName( "Proton" ),
    pionName(   "Pion" ),
  evSetup( &es ),
  prCollection( 0 ),
  piCollection( 0 ), 
   vCollection( v0Collection ),
   rCollection( 0 ),
  sList( searchList )  {
    ptSel = new BPHParticlePtSelect (  0.0 );
   etaSel = new BPHParticleEtaSelect( 10.0 );
  massSel = new BPHMassSelect( 0.0, 3.0 ) ;
  chi2Sel = 0;
  updated = false;
}


BPHLambda0ToPPiBuilder::BPHLambda0ToPPiBuilder(
               const edm::EventSetup& es,
               const std::vector<reco::VertexCompositePtrCandidate>*
                     vpCollection,
               const std::string& searchList ):
  protonName( "Proton" ),
    pionName(   "Pion" ),
  evSetup( &es ),
  prCollection( 0 ),
  piCollection( 0 ), 
   vCollection( 0 ),
   rCollection( vpCollection ),
  sList( searchList )  {
    ptSel = new BPHParticlePtSelect (  0.0 );
   etaSel = new BPHParticleEtaSelect( 10.0 );
  massSel = new BPHMassSelect( 0.0, 3.0 ) ;
  chi2Sel = 0;
  updated = false;
}

//--------------
// Destructor --
//--------------
BPHLambda0ToPPiBuilder::~BPHLambda0ToPPiBuilder() {
  delete   ptSel;
  delete  etaSel;
  delete massSel;
  delete chi2Sel;
}

//--------------
// Operations --
//--------------
vector<BPHPlusMinusConstCandPtr> BPHLambda0ToPPiBuilder::build(){

  if ( updated ) return lambda0List;
  
  if ( ( prCollection != 0 ) &&
       ( piCollection != 0 ) ) buildFromBPHGenericCollection();
  else
  if (    vCollection != 0   ) buildFromV0( vCollection );
  else
  if (    rCollection != 0   ) buildFromV0( rCollection );

  updated = true;
  return lambda0List;

}

void BPHLambda0ToPPiBuilder::buildFromBPHGenericCollection() {

  lambda0List.clear();

  BPHRecoBuilder bLambda0( *evSetup );
  bLambda0.add( protonName, prCollection, BPHParticleMasses::protonMass,
                                          BPHParticleMasses::protonMSigma );
  bLambda0.add(   pionName, piCollection, BPHParticleMasses::  pionMass,
                                          BPHParticleMasses::  pionMSigma );
  bLambda0.filter( protonName, *ptSel );
  bLambda0.filter(   pionName, *ptSel );
  bLambda0.filter( protonName, *etaSel );
  bLambda0.filter(   pionName, *etaSel );

  bLambda0.filter( *massSel );
  if ( chi2Sel != 0 )
  bLambda0.filter( *chi2Sel );

  lambda0List = BPHPlusMinusCandidate::build( bLambda0, protonName, pionName );

  return;

}

template <class T>
void BPHLambda0ToPPiBuilder::buildFromV0( const T* v0Collection ) {

  int iv0;
  int nv0 = v0Collection->size();
  lambda0List.clear();
  lambda0List.reserve( nv0 );
  
  // cycle over v0Collection
  for ( iv0 = 0; iv0 < nv0; ++iv0 ) {

    const typename T::value_type& lv0 = v0Collection->at( iv0 );

    // every reco::VertexCompositeCandidate must have exactly two daughters
    if ( lv0.numberOfDaughters() != 2 ) continue;

    // filters
    if( !ptSel ->accept( *lv0.daughter( 0 ) ) ) continue;
    if( !etaSel->accept( *lv0.daughter( 0 ) ) ) continue;
    if( !ptSel ->accept( *lv0.daughter( 1 ) ) ) continue;
    if( !etaSel->accept( *lv0.daughter( 1 ) ) ) continue;

    BPHPlusMinusCandidatePtr plX( new BPHPlusMinusCandidate( evSetup ) );
    BPHPlusMinusCandidatePtr plY( new BPHPlusMinusCandidate( evSetup ) );
    BPHPlusMinusCandidate* lambdaX = plX.get();
    BPHPlusMinusCandidate* lambdaY = plY.get();
    lambdaX->add( protonName, lv0.daughter( 0 ), sList,
                  BPHParticleMasses::protonMass );
    lambdaX->add(   pionName, lv0.daughter( 1 ), sList,
                  BPHParticleMasses::  pionMass );
    lambdaY->add( protonName, lv0.daughter( 1 ), sList,
                  BPHParticleMasses::protonMass );
    lambdaY->add(   pionName, lv0.daughter( 0 ), sList,
                  BPHParticleMasses::  pionMass );
    BPHPlusMinusCandidatePtr* pp0 = 0;

    // check which daughter has proton mass or choose the nearest V0 mass
    float mv0 = lv0.mass();
    if ( lv0.daughter( 0 )->mass() > 0.5 ) pp0 = &plX;
    else
    if ( lv0.daughter( 1 )->mass() > 0.5 ) pp0 = &plY;
    else pp0 = ( fabs( mv0 - lambdaX->mass() ) <
                 fabs( mv0 - lambdaY->mass() ) ? &plX : &plY );

    BPHPlusMinusCandidate* lambda0 = pp0->get();
    if ( !massSel->accept( *lambda0 ) ) continue;
    if ( ( chi2Sel != 0 ) &&
         ( !chi2Sel->accept( *lambda0 ) ) ) continue;

    const_cast<pat::CompositeCandidate&>( lambda0->composite() ).
                                          addUserFloat( "V0Mass", mv0 );
    lambda0List.push_back( *pp0 );

  }

  return;

}


/// set cuts
void BPHLambda0ToPPiBuilder::setPtMin( double pt ) {
  updated = false;
  ptSel->setPtMin( pt );
  return;
}


void BPHLambda0ToPPiBuilder::setEtaMax( double eta ) {
  updated = false;
  etaSel->setEtaMax( eta );
  return;
}


void BPHLambda0ToPPiBuilder::setMassMin( double m ) {
  updated = false;
  massSel->setMassMin( m );
  return;
}


void BPHLambda0ToPPiBuilder::setMassMax( double m ) {
  updated = false;
  massSel->setMassMax( m );
  return;
}


void BPHLambda0ToPPiBuilder::setProbMin( double p ) {
  updated = false;
  delete chi2Sel;
  chi2Sel = ( p < 0.0 ? 0 : new BPHChi2Select( p ) );
  return;
}


void BPHLambda0ToPPiBuilder::setConstr( double mass, double sigma ) {
  updated = false;
  cMass  = mass;
  cSigma = sigma;
  return;
}

/// get current cuts
double BPHLambda0ToPPiBuilder::getPtMin() const {
  return ptSel->getPtMin();
}


double BPHLambda0ToPPiBuilder::getEtaMax() const {
  return etaSel->getEtaMax();
}


double BPHLambda0ToPPiBuilder::getMassMin() const {
  return massSel->getMassMin();
}


double BPHLambda0ToPPiBuilder::getMassMax() const {
  return massSel->getMassMax();
}


double BPHLambda0ToPPiBuilder::getProbMin() const {
  return ( chi2Sel == 0 ? -1.0 : chi2Sel->getProbMin() );
}


double BPHLambda0ToPPiBuilder::getConstrMass() const {
  return cMass;
}


double BPHLambda0ToPPiBuilder::getConstrSigma() const {
  return cSigma;
}

