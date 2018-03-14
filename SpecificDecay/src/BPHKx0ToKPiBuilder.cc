/*
 *  See header file for a description of this class.
 *
 *  \author Paolo Ronchese INFN Padova
 *
 */

//-----------------------
// This Class' Header --
//-----------------------
#include "BPHAnalysis/SpecificDecay/interface/BPHKx0ToKPiBuilder.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BPHAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "BPHAnalysis/RecoDecay/interface/BPHPlusMinusCandidate.h"
#include "BPHAnalysis/RecoDecay/interface/BPHTrackReference.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticlePtSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticleEtaSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHMassSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHMassSymSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHChi2Select.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticleMasses.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"

/*
// old code left for example/reference
class BPHKx0ToKPiBuilder::DZSelect: public BPHRecoSelect {
 public:
  DZSelect( const std::string& nKaon, float dzMax,
            float mMin, float mMax ): name( nKaon ),
				      min( mMin * mMin ),
				      max( mMax * mMax ),
                                      dzm ( dzMax ) {}
  virtual ~DZSelect() {}
  virtual bool accept( const reco::Candidate& cand,
                       const BPHRecoBuilder* build ) const {
    const reco::Candidate& kaon = *get( name, build );
    const reco::Candidate::LorentzVector p4k = kaon.p4();
    const reco::Candidate::LorentzVector p4p = cand.p4();
    const float px = p4k.px() + p4p.px();
    const float py = p4k.py() + p4p.py();
    const float pz = p4k.pz() + p4p.pz();
    const float en = kaon.energy() + cand.energy();
    const float m2 = ( en * en ) - ( ( px * px ) + ( py * py ) + ( pz * pz ) );
    if ( m2 < min ) return false;
    if ( m2 > max ) return false;
    const reco::Track* kp = BPHTrackReference::getTrack( cand, "cfhp" );
    const reco::Track* kk = BPHTrackReference::getTrack( kaon, "cfhp" );
    if ( kp == 0 ) return false;
    if ( kk == 0 ) return false;
    return ( fabs( kk->dz() - kp->dz() ) < dzm );
  }
  void setMassMin( float m ) { min = m * m; }
  void setMassMax( float m ) { max = m * m; }
 private:
  std::string name;
  float min;
  float max;
  float dzm;
};
*/

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
BPHKx0ToKPiBuilder::BPHKx0ToKPiBuilder(
               const edm::EventSetup& es,
               const BPHRecoBuilder::BPHGenericCollection* pCollection,
               const BPHRecoBuilder::BPHGenericCollection* nCollection,
               BPHKx0ToKPiBuilder::buildMode mode ):
  kaonName( "Kaon" ),
  pionName( "Pion" ),
  evSetup( &es ),
  posCollection( pCollection ),
  negCollection( nCollection ) {
  ptMin = 0.7;
  etaMax = 10.0;
  massMin = 0.75;
  massMax = 1.05;
  probMin = 0.0;
// old code left for example/reference
//    ptSel = new BPHParticlePtSelect (  0.7 );
//   etaSel = new BPHParticleEtaSelect( 10.0 );
//  massSel = new BPHMassSelect( 0.75, 1.05 );
//  chi2Sel = new BPHChi2Select( 0.0 );
//    dzSel = new DZSelect( kaonName, 1.0, massSel->getMassMin(),
//                                         massSel->getMassMax() );
  updated = false;
}

//--------------
// Destructor --
//--------------
BPHKx0ToKPiBuilder::~BPHKx0ToKPiBuilder() {
// old code left for example/reference
//  delete   ptSel;
//  delete  etaSel;
//  delete massSel;
//  delete chi2Sel;
//  delete   dzSel;
}

//--------------
// Operations --
//--------------
vector<BPHPlusMinusConstCandPtr> BPHKx0ToKPiBuilder::build() {

  if ( updated ) return kx0List;

/*
// old code left for example/reference
  BPHRecoBuilder bKx0( *evSetup );
  bKx0.add( kaonName, pCollection, BPHParticleMasses::kaonMass,
                                   BPHParticleMasses::kaonMSigma );
  bKx0.add( pionName, nCollection, BPHParticleMasses::pionMass,
                                   BPHParticleMasses::pionMSigma );
  bKx0.filter( kaonName, * ptSel );
  bKx0.filter( pionName, * ptSel );
  bKx0.filter( kaonName, *etaSel );
  bKx0.filter( pionName, *etaSel );
  bKx0.filter( pionName, * dzSel );

  BPHMassSymSelect mTmpSel( kaonName, pionName, massSel );
  bKx0.filter( mTmpSel );

  vector<BPHPlusMinusConstCandPtr>
  tmpList = BPHPlusMinusCandidate::build( bKx0, kaonName, pionName );

  int ikx;
  int nkx = tmpList.size();
  kx0List.clear();
  kx0List.reserve( nkx );
  BPHPlusMinusConstCandPtr pxt( 0 );
  for ( ikx = 0; ikx < nkx; ++ikx ) {
    BPHPlusMinusConstCandPtr& px0 = tmpList[ikx];
    BPHPlusMinusCandidatePtr  pxb( new BPHPlusMinusCandidate( evSetup ) );
    const
    BPHPlusMinusCandidate* kx0 = px0.get();
    BPHPlusMinusCandidate* kxb = pxb.get();
    kxb->add( pionName, kx0->originalReco( kx0->getDaug( kaonName ) ),
              BPHParticleMasses::pionMass );
    kxb->add( kaonName, kx0->originalReco( kx0->getDaug( pionName ) ),
              BPHParticleMasses::kaonMass );
    if ( fabs( kx0->composite().mass() - BPHParticleMasses::kx0Mass ) <
         fabs( kxb->composite().mass() - BPHParticleMasses::kx0Mass ) )
         pxt = px0;
    else pxt = pxb;
    if ( !massSel->accept( *pxt ) ) continue;
    if ( !chi2Sel->accept( *pxt ) ) continue;
    kx0List.push_back( pxt );
  }
*/

  kx0List.clear();

  // extract basic informations from input collections

  class Particle {
   public:
    Particle( const reco::Candidate* c,
              const reco::Track* tk,
              double x,
              double y,
              double z,
              double p,
              double k ): cand( c ), track( tk ),
                          px( x ), py( y ), pz( z ), ePion( p ), eKaon( k ) {}
    const reco::Candidate* cand;
    const reco::Track* track;
    double px;
    double py;
    double pz;
    double ePion;
    double eKaon;
  };

  vector<Particle> pList;
  vector<Particle> nList;

  int nPos = posCollection->size();
  int nNeg = negCollection->size();

  pList.reserve( nPos );
  nList.reserve( nNeg );

  // filter input collections

  int iPos;
  int iNeg;

  for ( iPos = 0; iPos < nPos; ++iPos ) {
    const reco::Candidate& cand = posCollection->get( iPos );
    if ( cand.charge() != +1 ) continue;
    const reco::Candidate::LorentzVector p4 = cand.p4();
    if ( p4.pt () <  ptMin ) continue;
    if ( p4.eta() > etaMax ) continue;
    const reco::Track* tk = BPHTrackReference::getTrack( cand, "cfhp" );
    if ( tk == 0 ) continue;
    double px = p4.px();
    double py = p4.py();
    double pz = p4.pz();
    double p2 = ( px * px ) + ( py * py ) + ( pz * pz );
    pList.push_back( Particle( &cand, tk, px, py, pz,
                               sqrt( p2 + ( BPHParticleMasses::pionMass *
                                            BPHParticleMasses::pionMass ) ),
                               sqrt( p2 + ( BPHParticleMasses::kaonMass *
                                            BPHParticleMasses::kaonMass ) ) ) );
  }

  for ( iNeg = 0; iNeg < nNeg; ++iNeg ) {
    const reco::Candidate& cand = negCollection->get( iNeg );
    if ( cand.charge() != -1 ) continue;
    const reco::Candidate::LorentzVector p4 = cand.p4();
    if ( p4.pt () <  ptMin ) continue;
    if ( p4.eta() > etaMax ) continue;
    const reco::Track* tk = BPHTrackReference::getTrack( cand, "cfhp" );
    if ( tk == 0 ) continue;
    double px = p4.px();
    double py = p4.py();
    double pz = p4.pz();
    double p2 = ( px * px ) + ( py * py ) + ( pz * pz );
    nList.push_back( Particle( &cand, tk, px, py, pz,
                               sqrt( p2 + ( BPHParticleMasses::pionMass *
                                            BPHParticleMasses::pionMass ) ),
                               sqrt( p2 + ( BPHParticleMasses::kaonMass *
                                            BPHParticleMasses::kaonMass ) ) ) );
  }

  // filter basic candidates

  nPos = pList.size();
  nNeg = nList.size();

  double mMin = massMin * 0.9;
  double mMax = massMax * 1.1;

  for ( iPos = 0; iPos < nPos; ++iPos ) {
    Particle& pc = pList[iPos];
    const reco::Track* pt = pc.track;
    double px = pc.px;
    double py = pc.py;
    double pz = pc.pz;
    double pp = pc.ePion;
    double pk = pc.eKaon;
    for ( iNeg = 0; iNeg < nNeg; ++ iNeg ) {
      Particle& nc = nList[iNeg];
      const reco::Track* nt = nc.track;
      if ( fabs( nt->dz() - pt->dz() ) > 1.0 ) continue;
      double nx = nc.px;
      double ny = nc.py;
      double nz = nc.pz;
      double np = nc.ePion;
      double nk = nc.eKaon;
      const float tx = px + nx;
      const float ty = py + ny;
      const float tz = pz + nz;
      const float ta = pp + nk;
      const float tb = pk + np;
      float ma = ( ta * ta ) -
               ( ( tx * tx ) +
                 ( ty * ty ) +
                 ( tz * tz ) );
      float mb = ( tb * tb ) -
               ( ( tx * tx ) +
                 ( ty * ty ) +
                 ( tz * tz ) );
      if ( ma > 0.0 ) ma = sqrt( ma );
      else            ma = 0.0;
      if ( mb > 0.0 ) mb = sqrt( mb );
      else            mb = 0.0;
      double mass = -1;
      if ( fabs( ma - BPHParticleMasses::kx0Mass ) <
           fabs( mb - BPHParticleMasses::kx0Mass ) ) mass = ma;
      else                                           mass = mb;
      if ( mass < mMin ) continue;
      if ( mass > mMax ) continue;
      BPHPlusMinusCandidatePtr kx0( new BPHPlusMinusCandidate( evSetup ) );
      kx0->add( kaonName, pc.cand, BPHParticleMasses::kaonMass );
      kx0->add( pionName, nc.cand, BPHParticleMasses::pionMass );
      float massKx0 = kx0->composite().mass();
      BPHPlusMinusCandidatePtr kxb( new BPHPlusMinusCandidate( evSetup ) );
      kxb->add( kaonName, nc.cand, BPHParticleMasses::kaonMass );
      kxb->add( pionName, pc.cand, BPHParticleMasses::pionMass );
      float massKxb = kxb->composite().mass();
      BPHPlusMinusCandidatePtr kxt( 0 );
      if ( fabs( massKx0 - BPHParticleMasses::kx0Mass ) <
           fabs( massKxb - BPHParticleMasses::kx0Mass ) ) {
        mass = massKx0;
        kxt = kx0;
      }
      else {
        mass = massKxb;
        kxt = kxb;
      }
      if ( mass < massMin ) continue;
      if ( mass > massMax ) continue;
      const reco::Vertex& v = kxt->vertex();
      if (  v.isFake () ) continue;
      if ( !v.isValid() ) continue;
      if ( TMath::Prob( v.chi2(), lround( v.ndof() ) ) < probMin ) continue;
      kx0List.push_back( kxt );
    }
  }

  updated = true;
  return kx0List;

}

/// set cuts
void BPHKx0ToKPiBuilder::setPtMin( double pt ) {
  updated = false;
  ptMin = pt;
// old code left for example/reference
//  ptSel->setPtMin( pt );
  return;
}


void BPHKx0ToKPiBuilder::setEtaMax( double eta ) {
  updated = false;
  etaMax = eta;
// old code left for example/reference
//  etaSel->setEtaMax( eta );
  return;
}


void BPHKx0ToKPiBuilder::setMassMin( double m ) {
  updated = false;
  massMin = m;
// old code left for example/reference
//  massSel->setMassMin( m );
//    dzSel->setMassMin( m );
  return;
}


void BPHKx0ToKPiBuilder::setMassMax( double m ) {
  updated = false;
  massMax = m;
// old code left for example/reference
//  massSel->setMassMax( m );
//    dzSel->setMassMax( m );
  return;
}


void BPHKx0ToKPiBuilder::setProbMin( double p ) {
  updated = false;
  probMin = p;
// old code left for example/reference
//  chi2Sel->setProbMin( p );
  return;
}


// mass constraint: foreseen but actually not implemented
void BPHKx0ToKPiBuilder::setConstr( double mass, double sigma ) {
  updated = false;
  cMass  = mass;
  cSigma = sigma;
  return;
}

/// get current cuts
double BPHKx0ToKPiBuilder::getPtMin() const {
  return ptMin;
// old code left for example/reference
//  return ptSel->getPtMin();
}


double BPHKx0ToKPiBuilder::getEtaMax() const {
  return etaMax;
// old code left for example/reference
//  return etaSel->getEtaMax();
}


double BPHKx0ToKPiBuilder::getMassMin() const {
  return massMin;
// old code left for example/reference
//  return massSel->getMassMin();
}


double BPHKx0ToKPiBuilder::getMassMax() const {
  return massMax;
// old code left for example/reference
//  return massSel->getMassMax();
}


double BPHKx0ToKPiBuilder::getProbMin() const {
  return probMin;
// old code left for example/reference
//  return chi2Sel->getProbMin();
}


// mass constraint: foreseen but actually not implemented
double BPHKx0ToKPiBuilder::getConstrMass() const {
  return cMass;
}


// mass constraint: foreseen but actually not implemented
double BPHKx0ToKPiBuilder::getConstrSigma() const {
  return cSigma;
}

