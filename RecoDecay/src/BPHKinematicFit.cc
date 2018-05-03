/*
 *  See header file for a description of this class.
 *
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
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackMassKinematicConstraint.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

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
 oldKPs( true ),
 oldFit( true ),
 oldMom( true ),
 kinTree( 0 ) {
}


BPHKinematicFit::BPHKinematicFit( const BPHKinematicFit* ptr ):
 BPHDecayVertex( ptr, 0 ),
 massConst( -1.0 ),
 massSigma( -1.0 ),
 oldKPs( true ),
 oldFit( true ),
 oldMom( true ),
 kinTree( 0 ) {
  map<const reco::Candidate*,const reco::Candidate*> iMap;
  const vector<const reco::Candidate*>& daug = daughters();
  const vector<Component>& list = ptr->componentList();
  int i;
  int n = daug.size();
  for ( i = 0; i < n; ++i ) {
    const reco::Candidate* cand = daug[i];
    iMap[originalReco( cand )] = cand;
  }
  for ( i = 0; i < n; ++i ) {
    const Component& c = list[i];
    dMSig[iMap[c.cand]] = c.msig;
  }
  const vector<BPHRecoConstCandPtr>& dComp = daughComp();
  int j;
  int m = dComp.size();
  for ( j = 0; j < m; ++j ) {
    const BPHRecoCandidate* rc = dComp[j].get();
    const map<const reco:: Candidate*,double>& dMap = rc->dMSig;
    const map<const BPHRecoCandidate*,bool  >& cMap = rc->cKinP;
    dMSig.insert( dMap.begin(), dMap.end() );
    cKinP.insert( cMap.begin(), cMap.end() );
    cKinP[rc] = false;
  }
}

//--------------
// Destructor --
//--------------
BPHKinematicFit::~BPHKinematicFit() {
}

//--------------
// Operations --
//--------------
void BPHKinematicFit::setConstraint( double mass, double sigma ) {
  oldFit = oldMom = true;
  massConst = mass;
  massSigma = sigma;
  return;
}


double BPHKinematicFit::constrMass() const {
  return massConst;
}


double BPHKinematicFit::constrSigma() const {
  return massSigma;
}


void BPHKinematicFit::setIndependentFit( const string& name, bool flag ) {
  string::size_type pos = name.find( "/" );
  if ( pos != string::npos ) {
    edm::LogPrint( "WrongRequest" )
                << "BPHKinematicFit::setIndependentFit: "
                << "cascade decay specification not admitted " << name;
    return;
  }
  const BPHRecoCandidate* comp = getComp( name ).get();
  if ( comp == 0 ) {
    edm::LogPrint( "ParticleNotFound" )
                << "BPHKinematicFit::setIndependentFit: "
                << name << " daughter not found";
    return;
  }
  oldKPs = oldFit = oldMom = true;
  cKinP[comp] = flag;
  return;
}


const vector<RefCountedKinematicParticle>& BPHKinematicFit::kinParticles()
                                                            const {
  if ( oldKPs ) buildParticles();
  return allParticles;
}


vector<RefCountedKinematicParticle> BPHKinematicFit::kinParticles(
                                    const vector<string>& names ) const {
  if ( oldKPs ) buildParticles();
  vector<RefCountedKinematicParticle> plist;
  unsigned int m = allParticles.size();
  if ( m != numParticles() ) return plist;
  set<RefCountedKinematicParticle> pset;
  int i;
  int n = names.size();
  plist.reserve( m );
  for ( i = 0; i < n; ++i ) {
    const string& pname = names[i];
    if ( pname == "*" ) {
      int j = m;
      while ( j ) insertParticle( allParticles[--j], plist, pset );
      break;
    }
    string::size_type pos = pname.find( "/" );
    if ( pos != string::npos ) getParticles( pname.substr( 0, pos ),
                                             pname.substr( pos + 1 ),
                                             plist, pset );
    else                       getParticles( pname, "",
                                             plist, pset );
  }
  return plist;
}


const RefCountedKinematicTree& BPHKinematicFit::kinematicTree() const {
  if ( oldFit ) return kinematicTree( "", massConst, massSigma );
  return kinTree;
}


const RefCountedKinematicTree& BPHKinematicFit::kinematicTree(
                               const string& name,
                               double mass, double sigma ) const {
  if ( mass  < 0 ) return kinematicTree( name );
  if ( sigma < 0 ) return kinematicTree( name, mass );
  ParticleMass mc = mass;
  MassKinematicConstraint kinConst( mc, sigma );
  return kinematicTree( name, &kinConst );
}


const RefCountedKinematicTree& BPHKinematicFit::kinematicTree(
                               const string& name,
                               double mass ) const {
  if ( mass < 0 ) {
    kinTree = RefCountedKinematicTree( 0 );
    oldFit = false;
    return kinTree;
  }
  int nn = daughFull().size();
  ParticleMass mc = mass;
  if ( nn == 2 ) {
    TwoTrackMassKinematicConstraint   kinConst( mc );
    return kinematicTree( name, &kinConst );
  }
  else {
    MultiTrackMassKinematicConstraint kinConst( mc, nn );
    return kinematicTree( name, &kinConst );
  }
}


const RefCountedKinematicTree& BPHKinematicFit::kinematicTree(
                               const string& name ) const {
  KinematicConstraint* kc = 0;
  return kinematicTree( name, kc );
}


const RefCountedKinematicTree& BPHKinematicFit::kinematicTree(
                               const string& name,
                               KinematicConstraint* kc ) const {
  kinTree = RefCountedKinematicTree( 0 );
  oldFit = false;
  kinParticles();
  if ( allParticles.size() != numParticles() ) return kinTree;
  vector<RefCountedKinematicParticle> kComp;
  vector<RefCountedKinematicParticle> kTail;
  splitKP( name, &kComp, &kTail );
  if ( kComp.size() == 0 ) {
    edm::LogPrint( "ParticleNotFound" )
                << "BPHKinematicFit::kinematicTree: "
                << name << " daughter not found";
    return kinTree;
  }
  try {
    KinematicParticleVertexFitter vtxFitter;
    RefCountedKinematicTree compTree = vtxFitter.fit( kComp );
    if ( compTree->isEmpty() ) return kinTree;
    if ( kc != 0 ) {
      KinematicParticleVertexFitter vtxFitter;
      KinematicParticleFitter kinFitter;
      compTree = kinFitter.fit( kc, compTree );
      if ( compTree->isEmpty() ) return kinTree;
    }
    compTree->movePointerToTheTop();
    if ( kTail.size() ) {
      RefCountedKinematicParticle compPart = compTree->currentParticle();
      if ( !compPart->currentState().isValid() ) return kinTree;
      kTail.push_back( compPart );
      kinTree = vtxFitter.fit( kTail );
    }
    else {
      kinTree = compTree;
    }
  }
  catch ( std::exception e ) {
    edm::LogPrint( "FitFailed" )
                << "BPHKinematicFit::kinematicTree: "
                << "kin fit reset";
    kinTree = RefCountedKinematicTree( 0 );
  }
  return kinTree;
}


const RefCountedKinematicTree& BPHKinematicFit::kinematicTree(
                               const string& name,
                               MultiTrackKinematicConstraint* kc ) const {
  kinTree = RefCountedKinematicTree( 0 );
  oldFit = false;
  kinParticles();
  if ( allParticles.size() != numParticles() ) return kinTree;
  vector<RefCountedKinematicParticle> kPart;
  splitKP( name, &kPart );
  if ( kPart.size() == 0 ) {
    edm::LogPrint( "ParticleNotFound" )
                << "BPHKinematicFit::kinematicTree: "
                << name << " daughter not found";
    return kinTree;
  }
  try {
    KinematicConstrainedVertexFitter cvf;
    kinTree = cvf.fit( kPart, kc );
  }
  catch ( std::exception e ) {
    edm::LogPrint( "FitFailed" )
                << "BPHKinematicFit::kinematicTree: "
                << "kin fit reset";
    kinTree = RefCountedKinematicTree( 0 );
  }
  return kinTree;
}


void BPHKinematicFit::resetKinematicFit() const {
  oldKPs = oldFit = oldMom = true;
  return;
}


bool BPHKinematicFit::isEmpty() const {
  kinematicTree();
  if ( kinTree.get() == 0 ) return true;
  return kinTree->isEmpty();
}


bool BPHKinematicFit::isValidFit() const {
  const RefCountedKinematicParticle kPart = topParticle();
  if ( kPart.get() == 0 ) return false;
  return kPart->currentState().isValid();
}


const RefCountedKinematicParticle BPHKinematicFit::currentParticle() const {
  if ( isEmpty() ) return RefCountedKinematicParticle( 0 );
  return kinTree->currentParticle();
}


const RefCountedKinematicVertex BPHKinematicFit::currentDecayVertex() const {
  if ( isEmpty() ) return RefCountedKinematicVertex( 0 );
  return kinTree->currentDecayVertex();
}


const RefCountedKinematicParticle BPHKinematicFit::topParticle() const {
  if ( isEmpty() ) return RefCountedKinematicParticle( 0 );
  return kinTree->topParticle();
}


const RefCountedKinematicVertex BPHKinematicFit::topDecayVertex() const {
  if ( isEmpty() ) return RefCountedKinematicVertex( 0 );
  kinTree->movePointerToTheTop();
  return kinTree->currentDecayVertex();
}


ParticleMass BPHKinematicFit::mass() const {
  const RefCountedKinematicParticle kPart = topParticle();
  if ( kPart.get() == 0 ) return -1.0;
  const KinematicState kStat = kPart->currentState();
  if ( kStat.isValid() ) return kStat.mass();
  return -1.0;
}


const math::XYZTLorentzVector& BPHKinematicFit::p4() const {
  if ( oldMom ) fitMomentum();
  return totalMomentum;
}


void BPHKinematicFit::addK( const string& name,
                            const reco::Candidate* daug, 
                            double mass, double sigma ) {
  addK( name, daug, "cfhpmig", mass, sigma );
  return;
}


void BPHKinematicFit::addK( const string& name,
                            const reco::Candidate* daug, 
                            const string& searchList,
                            double mass, double sigma ) {
  addV( name, daug, searchList, mass );
  dMSig[daughters().back()] = sigma;
  return;
}


void BPHKinematicFit::addK( const string& name,
                            const BPHRecoConstCandPtr& comp ) {
  addV( name, comp );
  const map<const reco::Candidate*,double>& dMap = comp->dMSig;
  const map<const BPHRecoCandidate*,bool >& cMap = comp->cKinP;
  dMSig.insert( dMap.begin(), dMap.end() );
  cKinP.insert( cMap.begin(), cMap.end() );
  cKinP[comp.get()] = false;
  return;
}


void BPHKinematicFit::setNotUpdated() const {
  BPHDecayVertex::setNotUpdated();
  resetKinematicFit();
  return;
}


void BPHKinematicFit::buildParticles() const {
  kinMap.clear();
  kCDMap.clear();
  allParticles.clear();
  allParticles.reserve( daughFull().size() );
  addParticles( allParticles, kinMap, kCDMap );
  oldKPs = false;
  return;
}


void BPHKinematicFit::addParticles(
                      vector<RefCountedKinematicParticle>& kl,
                      map   <const reco::Candidate*,
                             RefCountedKinematicParticle>& km,
                      map   <const BPHRecoCandidate*,
                             RefCountedKinematicParticle>& cm ) const {
  const vector<const reco::Candidate*>& daug = daughters();
  KinematicParticleFactoryFromTransientTrack pFactory;
  int n = daug.size();
  float chi = 0.0;
  float ndf = 0.0;
  while ( n-- ) {
    const reco::Candidate* cand = daug[n];
    ParticleMass mass = cand->mass();
    float sigma = dMSig.find( cand )->second;
    if ( sigma < 0 ) sigma = 1.0e-7;
    reco::TransientTrack* tt = getTransientTrack( cand );
    if ( tt != 0 ) kl.push_back( km[cand] = pFactory.particle( *tt, 
                                            mass, chi, ndf, sigma ) );
  }
  const vector<BPHRecoConstCandPtr>& comp = daughComp();
  int m = comp.size();
  while ( m-- ) {
    const BPHRecoCandidate* cptr = comp[m].get();
    if ( cKinP.at( cptr ) ) {
      VirtualKinematicParticleFactory vFactory;
      if ( cptr->isEmpty() ) return;
      if ( !cptr->isValidFit() ) return;
      const RefCountedKinematicParticle part = cptr->topParticle();
      const RefCountedKinematicVertex   pvtx = cptr->topDecayVertex();
      float chi = pvtx->chiSquared();
      float dof = pvtx->degreesOfFreedom();
      kl.push_back( cm[cptr] = vFactory.particle( part->currentState(),
                                                  chi, dof, part ) );
    }
    else {
      cptr->addParticles( kl, km, cm );
    }
  }
  return;
}


void BPHKinematicFit::getParticles( const string& moth, const string& daug,
                      vector<RefCountedKinematicParticle>& kl,
                      set   <RefCountedKinematicParticle>& ks ) const {
  const BPHRecoCandidate* cptr = getComp( moth ).get();
  if ( cptr != 0 ) {
    if ( cKinP.at( cptr ) ) {
      insertParticle( kCDMap[cptr], kl, ks );
    }
    else {
      vector<string> list;
      if ( daug != "" ) {
        list.push_back( daug );
      }
      else {
        const vector<string>& dNames = cptr->daugNames();
        const vector<string>& cNames = cptr->compNames();
        list.insert( list.end(), dNames.begin(), dNames.end() );
        list.insert( list.end(), cNames.begin(), cNames.end() );
      }
      getParticles( moth, list, kl, ks );
    }
    return;
  }
  const reco::Candidate* dptr = getDaug( moth );
  if ( dptr != 0 ) {
    insertParticle( kinMap[dptr], kl, ks );
    return;
  }
  edm::LogPrint( "ParticleNotFound" )
              << "BPHKinematicFit::getParticles: "
              << moth << " not found";
  return;
}


void BPHKinematicFit::getParticles( const string& moth, const
                      vector<string>& daug,
                      vector<RefCountedKinematicParticle>& kl,
                      set   <RefCountedKinematicParticle>& ks ) const {
  int i;
  int n = daug.size();
  for ( i = 0; i < n; ++i ) {
    const string& name = daug[i];
    string::size_type pos = name.find( "/" );
    if ( pos != string::npos )
         getParticles( moth + "/" + name.substr( 0, pos ), 
                       name.substr( pos + 1 ),
                       kl, ks );
    else getParticles( moth + "/" + name, "",
                       kl, ks );
  }
  return;
}


unsigned int BPHKinematicFit::numParticles(
                              const BPHKinematicFit* cand ) const {
  if ( cand == 0 ) cand = this;
  unsigned int n = cand->daughters().size();
  const vector<string>& cnames = cand->compNames();
  int i = cnames.size();
  while ( i ) {
    const BPHRecoCandidate* comp = cand->getComp( cnames[--i] ).get();
    if ( cKinP.at( comp ) ) ++n;
    else                    n += numParticles( comp );
  }
  return n;
}


void BPHKinematicFit::insertParticle( RefCountedKinematicParticle& kp,
                      vector<RefCountedKinematicParticle>& kl,
                      set   <RefCountedKinematicParticle>& ks ) {
  if ( ks.find( kp ) != ks.end() ) return;
  kl.push_back( kp );
  ks.insert   ( kp );
  return;
}


void BPHKinematicFit::splitKP( const string& name,
                      vector<RefCountedKinematicParticle>* kComp,
                      vector<RefCountedKinematicParticle>* kTail ) const {
  kComp->clear();
  if ( kTail == 0 ) kTail = kComp;
  else              kTail->clear();
  kinParticles();
  if ( name == "" ) {
    *kComp = allParticles;
    return;
  }
  const reco:: Candidate* daug = getDaug( name );
  const BPHRecoCandidate* comp = getComp( name ).get();
  int ns;
  if ( daug != 0 ) {
    ns = 1;
  }
  else
  if ( comp != 0 ) {
    ns = ( cKinP.at( comp ) ? 1 : numParticles( comp ) );
  }
  else {
    edm::LogPrint( "ParticleNotFound" )
                << "BPHKinematicFit::splitKP: "
                << name << " daughter not found";
    *kTail = allParticles;
    return;
  }
  vector<string> nfull( 2 );
  nfull[0] = name;
  nfull[1] = "*";
  vector<RefCountedKinematicParticle> kPart = kinParticles( nfull );
  vector<RefCountedKinematicParticle>::const_iterator iter = kPart.begin();
  vector<RefCountedKinematicParticle>::const_iterator imid = iter + ns;
  vector<RefCountedKinematicParticle>::const_iterator iend = kPart.end();
  kComp->insert( kComp->end(), iter, imid );
  kTail->insert( kTail->end(), imid, iend );
  return;
}


void BPHKinematicFit::fitMomentum() const {
  if ( isValidFit() ) {
    const KinematicState& ks = currentParticle()->currentState();
    GlobalVector tm = ks.globalMomentum();
    double x = tm.x();
    double y = tm.y();
    double z = tm.z();
    double m = ks.mass();
    double e = sqrt( ( x * x ) + ( y * y ) + ( z * z ) + ( m * m ) );
    totalMomentum.SetPxPyPzE( x, y, z, e );
  }
  else {
    edm::LogPrint( "FitNotFound" )
                << "BPHKinematicFit::fitMomentum: "
                << "simple momentum sum computed";
    math::XYZTLorentzVector tm;
    const vector<const reco::Candidate*>& daug = daughters();
    int n = daug.size();
    while ( n-- ) tm += daug[n]->p4();
    const vector<BPHRecoConstCandPtr>& comp = daughComp();
    int m = comp.size();
    while ( m-- ) tm += comp[m]->p4();
    totalMomentum = tm;
  }
  oldMom = false;
  return;
}

