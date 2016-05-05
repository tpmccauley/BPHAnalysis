/*
 *  See header file for a description of this class.
 *
 *  $Date: 2015-07-03 17:02:36 $
 *  $Revision: 1.1 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//-----------------------
// This Class' Header --
//-----------------------
#include "BPHAnalysis/RecoDecay/interface/BPHKinematicFit.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BPHAnalysis/RecoDecay/interface/BPHRecoCandidate.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"

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
BPHKinematicFit::BPHKinematicFit():
  BPHDecayVertex( 0 ),
  massConst( -1.0 ),
  massSigma( -1.0 ),
  updatedFit( false ),
  updatedMom( false ),
  kinTree( 0 ) {
}


BPHKinematicFit::BPHKinematicFit( const BPHKinematicFit* ptr ):
 BPHDecayVertex( ptr, 0 ),
 massConst( -1.0 ),
 massSigma( -1.0 ),
 updatedFit( false ),
 updatedMom( false ),
 kinTree( 0 ) {
  std::map<const reco::Candidate*,
           const reco::Candidate*> iMap;
  const vector<const reco::Candidate*>& daug = daughters();
  const std::vector<Component>& list = ptr->BPHDecayMomentum::componentList();
  int i;
  int n = list.size();
  for ( i = 0; i < n; ++i ) {
    const reco::Candidate* cand = daug[i];
    iMap[originalReco( cand )] = cand;
  }
  for ( i = 0; i < n; ++i ) {
    const Component& c = list[i];
    dMSig[iMap[c.cand]] = c.msig;
  }
}

//--------------
// Destructor --
//--------------
BPHKinematicFit::~BPHKinematicFit() {
  delete kinTree;
}

//--------------
// Operations --
//--------------
void BPHKinematicFit::add( const std::string& name,
                           const reco::Candidate* daug, 
                           double mass, double sigma ) {
  add( name, daug, "cfhpmig", mass, sigma );
  return;
}


void BPHKinematicFit::add( const std::string& name,
                           const reco::Candidate* daug, 
                           const std::string& searchList,
                           double mass, double sigma ) {
  BPHDecayVertex::add( name, daug, searchList, mass );
  dMSig[daughters().back()] = sigma;
  return;
}


void BPHKinematicFit::setConstraint( double mass, double sigma ) {
  massConst = mass;
  massSigma = sigma;
  updatedFit = updatedMom = false;
  return;
}


double BPHKinematicFit::constrMass() const {
  return massConst;
}


double BPHKinematicFit::constrSigma() const {
  return massSigma;
}


const RefCountedKinematicTree& BPHKinematicFit::kinematicTree() const {
  if ( !updatedFit ) kinFit();
  return *kinTree;
}


const math::XYZTLorentzVector& BPHKinematicFit::p4() const {
  if ( !updatedMom ) fitMomentum();
  return totalMomentum;
}


void BPHKinematicFit::setNotUpdated() const {
  BPHDecayVertex::setNotUpdated();
  updatedFit = updatedMom = false;
  return;
}


void BPHKinematicFit::kinFit() const {
  delete kinTree;
  kinTree = new RefCountedKinematicTree;
  if ( massConst < 0 ) return;
  const std::vector<const reco::Candidate*>& daug = daughFull();
  KinematicParticleFactoryFromTransientTrack pFactory;
  std::vector<RefCountedKinematicParticle> allParticles;
  int n = daug.size();
  float chi = 0.0;
  float ndf = 0.0;
  while ( n-- ) {
    const reco::Candidate* cand = daug[n];
    ParticleMass m = cand->mass();
    float s = dMSig.find( cand )->second;
    if ( s < 0 ) s = 1.0e-10;
    reco::TransientTrack* tt = getTransientTrack( cand );
    if ( tt != 0 ) allParticles.push_back( pFactory.particle( *tt, 
                                           m, chi, ndf, s ) );
  }
  KinematicParticleVertexFitter vtxFitter;
  updatedFit = true;
  *kinTree = vtxFitter.fit( allParticles );
  RefCountedKinematicTree& kt = *kinTree;
  if ( kt->isEmpty() ) return;
  KinematicParticleFitter kinFitter;
  double mSig = ( massSigma < 0 ?  1.0e-10 : massSigma );
  MassKinematicConstraint kinConst( massConst, mSig );
  kt = kinFitter.fit( &kinConst, kt );
  if ( kt->isEmpty() ) return;
  kt->movePointerToTheTop();
  return;
}


void BPHKinematicFit::fitMomentum() const {
  if ( massConst < 0 ) {
    math::XYZTLorentzVector tm;
    const std::vector<const reco::Candidate*>& daug = daughters();
    int n = daug.size();
    while ( n-- ) tm += daug[n]->p4();
    const std::vector<BPHRecoConstCandPtr>& comp = daughComp();
    int m = comp.size();
    while ( m-- ) tm += comp[m]->p4();
    totalMomentum = tm;
  }
  else {
    const RefCountedKinematicTree& kt = kinematicTree();
    const KinematicState& ks = kt->currentParticle()->currentState();
    GlobalVector tm = ks.globalMomentum();
    double x = tm.x();
    double y = tm.y();
    double z = tm.z();
    double m = ks.mass();
    double e = sqrt( ( x * x ) + ( y * y ) + ( z * z ) + ( m * m ) );
    totalMomentum.SetPxPyPzE( x, y, z, e );
  }
  updatedMom = true;
  return;
}

