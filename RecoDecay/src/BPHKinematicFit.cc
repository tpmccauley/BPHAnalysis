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
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
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
 updatedKPs( false ),
 updatedFit( false ),
 updatedMom( false ),
 kinTree( 0 ) {
}


BPHKinematicFit::BPHKinematicFit( const BPHKinematicFit* ptr ):
 BPHDecayVertex( ptr, 0 ),
 massConst( -1.0 ),
 massSigma( -1.0 ),
 updatedKPs( false ),
 updatedFit( false ),
 updatedMom( false ),
 kinTree( 0 ) {
  map<const reco::Candidate*,const reco::Candidate*> iMap;
  const vector<const reco::Candidate*>& daug = daughters();
  const vector<Component>& list = ptr->BPHDecayMomentum::componentList();
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
void BPHKinematicFit::add( const string& name,
                           const reco::Candidate* daug, 
                           double mass, double sigma ) {
  add( name, daug, "cfhpmig", mass, sigma );
  return;
}


void BPHKinematicFit::add( const string& name,
                           const reco::Candidate* daug, 
                           const string& searchList,
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


const vector<RefCountedKinematicParticle>& BPHKinematicFit::kinParticles()
                                                            const {
  if ( !updatedKPs ) buildParticles();
  return allParticles;
}


vector<RefCountedKinematicParticle> BPHKinematicFit::kinParticles(
                                    const vector<string>& names ) const {
  if ( !updatedKPs ) buildParticles();
  set<RefCountedKinematicParticle> pset;
  vector<RefCountedKinematicParticle> plist;
  const vector<const reco::Candidate*>& daugs = daughFull();
  int i;
  int n = names.size();
  int m = daugs.size();
  plist.reserve( m );
  for ( i = 0; i < n; ++i ) {
    const string& pname = names[i];
    if ( pname == "*" ) {
      int j = m;
      while ( j-- ) {
        RefCountedKinematicParticle& kp = allParticles[j];
        if ( pset.find( kp ) != pset.end() ) continue;
        plist.push_back( kp );
        pset .insert   ( kp );
      }
      break;
    }
    map<const reco::Candidate*,
        RefCountedKinematicParticle>::const_iterator iter = kinMap.find(
                                                            getDaug( pname ) );
    map<const reco::Candidate*,
        RefCountedKinematicParticle>::const_iterator iend = kinMap.end();
    if ( iter != iend ) {
      const RefCountedKinematicParticle& kp = iter->second;
      if ( pset.find( kp ) != pset.end() ) continue;
      plist.push_back( kp );
      pset .insert   ( kp );
    }
    else {
      cout << "BPHKinematicFit::kinParticles: " << pname << " not found"
           << endl;
    }
  }
  return plist;
}


const RefCountedKinematicTree& BPHKinematicFit::kinematicTree() const {
  if ( !updatedFit ) kinFit();
  return *kinTree;
}


const RefCountedKinematicTree& BPHKinematicFit::kinematicTree(
                               const string& name,
                               MultiTrackKinematicConstraint* kc ) const {
  if ( !updatedFit ) kinFit( name, kc );
  return *kinTree;
}


void BPHKinematicFit::setKinematicFit( const RefCountedKinematicTree& kt ) {
  delete kinTree;
  kinTree = new RefCountedKinematicTree( kt );
  updatedFit = true;
  return;
}


void BPHKinematicFit::resetKinematicFit() {
  setNotUpdated();
  return;
}


const math::XYZTLorentzVector& BPHKinematicFit::p4() const {
  if ( !updatedMom ) fitMomentum();
  return totalMomentum;
}


void BPHKinematicFit::setNotUpdated() const {
  BPHDecayVertex::setNotUpdated();
  updatedKPs = updatedFit = updatedMom = false;
  return;
}


void BPHKinematicFit::buildParticles() const {
  kinMap.clear();
  allParticles.clear();
  const vector<const reco::Candidate*>& daug = daughFull();
  KinematicParticleFactoryFromTransientTrack pFactory;
  int n = daug.size();
  allParticles.reserve( n );
  float chi = 0.0;
  float ndf = 0.0;
  while ( n-- ) {
    const reco::Candidate* cand = daug[n];
    ParticleMass m = cand->mass();
    float s = dMSig.find( cand )->second;
    if ( s < 0 ) s = 1.0e-10;
    reco::TransientTrack* tt = getTransientTrack( cand );
    if ( tt != 0 ) allParticles.push_back( kinMap[cand] =
                                           pFactory.particle( *tt, 
                                           m, chi, ndf, s ) );
  }
  updatedKPs = true;
  return;
}


void BPHKinematicFit::kinFit() const {
  delete kinTree;
  kinTree = new RefCountedKinematicTree( 0 );
  updatedFit = true;
  if ( massConst < 0 ) return;
  kinParticles();
  if ( allParticles.size() != daughFull().size() ) return;
  try {
    KinematicParticleVertexFitter vtxFitter;
    *kinTree = vtxFitter.fit( allParticles );
    RefCountedKinematicTree& kt = *kinTree;
    if ( kt->isEmpty() ) return;
    KinematicParticleFitter kinFitter;
    double mSig = ( massSigma < 0 ?  1.0e-10 : massSigma );
    MassKinematicConstraint kinConst( massConst, mSig );
    kt = kinFitter.fit( &kinConst, kt );
    if ( kt->isEmpty() ) return;
    kt->movePointerToTheTop();
  }
  catch ( std::exception e ) {
    cout << "kin fit failed, reset" << endl;
    delete kinTree;
    kinTree = new RefCountedKinematicTree( 0 );
  }
  return;
}


void BPHKinematicFit::kinFit( const string& name,
                              MultiTrackKinematicConstraint* kc ) const {
  const BPHRecoCandidate* comp = getComp( name ).get();
  delete kinTree;
  kinTree = new RefCountedKinematicTree( 0 );
  updatedFit = true;
  kinParticles();
  if ( allParticles.size() != daughFull().size() ) return;
  const vector<string>& names = comp->daugNames();
  int nn = names.size();
  vector<string> nfull( nn + 1 );
  nfull[nn] = "*";
  while ( nn-- ) nfull[nn] = name + "/" + names[nn];
  try {
    KinematicConstrainedVertexFitter cvf;
    *kinTree = cvf.fit( kinParticles( nfull ), kc );
  }
  catch ( std::exception e ) {
    cout << "kin fit failed, reset" << endl;
    delete kinTree;
    kinTree = new RefCountedKinematicTree( 0 );
  }
  return;
}


void BPHKinematicFit::fitMomentum() const {
  if ( kinTree == 0 ) kinFit();
  const KinematicTree* kt = kinTree->get();
  if ( ( kt == 0 ) || ( kt->isEmpty() ) ) {
    if ( kt != 0 ) cout << "kin fit failed, simple momentum sum computed"
                        << endl;
    math::XYZTLorentzVector tm;
    const vector<const reco::Candidate*>& daug = daughters();
    int n = daug.size();
    while ( n-- ) tm += daug[n]->p4();
    const vector<BPHRecoConstCandPtr>& comp = daughComp();
    int m = comp.size();
    while ( m-- ) tm += comp[m]->p4();
    totalMomentum = tm;
  }
  else {
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

