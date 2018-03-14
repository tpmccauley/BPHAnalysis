/*
 *  See header file for a description of this class.
 *
 *  \author Paolo Ronchese INFN Padova
 *
 */

//-----------------------
// This Class' Header --
//-----------------------
#include "BPHAnalysis/SpecificDecay/interface/BPHPhiToKKBuilder.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BPHAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "BPHAnalysis/RecoDecay/interface/BPHPlusMinusCandidate.h"
#include "BPHAnalysis/RecoDecay/interface/BPHTrackReference.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticlePtSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticleEtaSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHMassSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHChi2Select.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticleMasses.h"

/*
// old code left for example/reference
class BPHPhiToKKBuilder::DZSelect: public BPHRecoSelect {
 public:
  DZSelect( const std::string& nKPos, float dzMax,
            float mMin, float mMax ): name( nKPos ),
				      min( mMin * mMin ),
				      max( mMax * mMax ),
                                      dzm ( dzMax ) {}
  virtual ~DZSelect() {}
  virtual bool accept( const reco::Candidate& cand,
                       const BPHRecoBuilder* build ) const {
    const reco::Candidate& cpos = *get( name, build );
    const reco::Candidate::LorentzVector p4p = cpos.p4();
    const reco::Candidate::LorentzVector p4n = cand.p4();
    const float px = p4p.px() + p4n.px();
    const float py = p4p.py() + p4n.py();
    const float pz = p4p.pz() + p4n.pz();
    const float en = cpos.energy() + cand.energy();
    const float m2 = ( en * en ) - ( ( px * px ) + ( py * py ) + ( pz * pz ) );
    if ( m2 < min ) return false;
    if ( m2 > max ) return false;
    const reco::Track* kn = BPHTrackReference::getTrack( cand, "cfhp" );
    const reco::Track* kp = BPHTrackReference::getTrack( cpos, "cfhp" );
    if ( kn == 0 ) return false;
    if ( kp == 0 ) return false;
    return ( fabs( kp->dz() - kn->dz() ) < dzm );
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
using namespace std;

//-------------------
// Initializations --
//-------------------


//----------------
// Constructors --
//----------------
BPHPhiToKKBuilder::BPHPhiToKKBuilder(
               const edm::EventSetup& es,
               const BPHRecoBuilder::BPHGenericCollection* pCollection,
               const BPHRecoBuilder::BPHGenericCollection* nCollection ):
  kPosName( "KPos" ),
  kNegName( "KNeg" ),
  evSetup( &es ),
  posCollection( pCollection ),
  negCollection( nCollection ) {
  ptMin = 0.7;
  etaMax = 10.0;
  massMin = 1.0;
  massMax = 1.04;
  probMin = 0.0;
// old code left for example/reference
//    ptSel = new BPHParticlePtSelect (  0.7 );
//   etaSel = new BPHParticleEtaSelect( 10.0 );
//  massSel = new BPHMassSelect( 1.0, 1.04 );
//  chi2Sel = new BPHChi2Select( 0.0 );
//    dzSel = new DZSelect( kPosName, 1.0, massSel->getMassMin(),
//                                         massSel->getMassMax() );
  updated = false;
}

//--------------
// Destructor --
//--------------
BPHPhiToKKBuilder::~BPHPhiToKKBuilder() {
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
vector<BPHPlusMinusConstCandPtr> BPHPhiToKKBuilder::build() {

  if ( updated ) return phiList;

/*
// old code left for example/reference
  BPHRecoBuilder bPhi( *evSetup );
  bPhi.add( kPosName, posCollection, BPHParticleMasses::kaonMass,
                                     BPHParticleMasses::kaonMSigma );
  bPhi.add( kNegName, negCollection, BPHParticleMasses::kaonMass,
                                     BPHParticleMasses::kaonMSigma );
  bPhi.filter( kPosName, * ptSel );
  bPhi.filter( kNegName, * ptSel );
  bPhi.filter( kPosName, *etaSel );
  bPhi.filter( kNegName, *etaSel );
  bPhi.filter( kNegName, * dzSel );

  bPhi.filter( *massSel );
  bPhi.filter( *chi2Sel );

  phiList = BPHPlusMinusCandidate::build( bPhi, kPosName, kNegName );
*/

  phiList.clear();

  // extract basic informations from input collections

  class Particle {
   public:
    Particle( const reco::Candidate* c,
          const reco::Track* tk,
          double x,
          double y,
          double z,
          double e ): cand( c ), track( tk ),
                      px( x ), py( y ), pz( z ), energy( e ) {}
    const reco::Candidate* cand;
    const reco::Track* track;
    double px;
    double py;
    double pz;
    double energy;
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
    pList.push_back( Particle( &cand, tk, px, py, pz,
                               sqrt( ( px * px ) + ( py * py ) + ( pz * pz ) +
                                     ( BPHParticleMasses::kaonMass *
                                       BPHParticleMasses::kaonMass ) ) ) );
  }

  for ( iNeg = 0; iNeg < nNeg; ++ iNeg ) {
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
    nList.push_back( Particle( &cand, tk, px, py, pz,
                               sqrt( ( px * px ) + ( py * py ) + ( pz * pz ) +
                                     ( BPHParticleMasses::kaonMass *
                                       BPHParticleMasses::kaonMass ) ) ) );
  }

  nPos = pList.size();
  nNeg = nList.size();

  double mSqMin = massMin * massMin * 0.9;
  double mSqMax = massMax * massMax * 1.1;

  for ( iPos = 0; iPos < nPos; ++iPos ) {
    Particle& pc = pList[iPos];
    const reco::Track* pt = pc.track;
    double px = pc.px;
    double py = pc.py;
    double pz = pc.pz;
    double pe = pc.energy;
    for ( iNeg = 0; iNeg < nNeg; ++ iNeg ) {
      Particle& nc = nList[iNeg];
      const reco::Track* nt = nc.track;
      if ( fabs( nt->dz() - pt->dz() ) > 1.0 ) continue;
      double nx = nc.px;
      double ny = nc.py;
      double nz = nc.pz;
      double ne = nc.energy;
      const float tx = px + nx;
      const float ty = py + ny;
      const float tz = pz + nz;
      const float te = pe + ne;
      const float m2 = ( te * te ) -
                     ( ( tx * tx ) +
                       ( ty * ty ) +
                       ( tz * tz ) );
      if ( m2 < mSqMin ) continue;
      if ( m2 > mSqMax ) continue;
      BPHPlusMinusCandidatePtr phi( new BPHPlusMinusCandidate( evSetup ) );
      phi->add( kPosName, pc.cand, BPHParticleMasses::kaonMass );
      phi->add( kNegName, nc.cand, BPHParticleMasses::kaonMass );
      double mass = phi->composite().mass();
      if ( mass < massMin ) continue;
      if ( mass > massMax ) continue;
      const reco::Vertex& v = phi->vertex();
      if (  v.isFake () ) continue;
      if ( !v.isValid() ) continue;
      if ( TMath::Prob( v.chi2(), lround( v.ndof() ) ) < probMin ) continue;
      phiList.push_back( phi );
    }
  }

  updated = true;
  return phiList;

}

/// set cuts
void BPHPhiToKKBuilder::setPtMin( double pt ) {
  updated = false;
  ptMin = pt;
//  ptSel->setPtMin( pt );
  return;
}


void BPHPhiToKKBuilder::setEtaMax( double eta ) {
  updated = false;
  etaMax = eta;
//  etaSel->setEtaMax( eta );
  return;
}


void BPHPhiToKKBuilder::setMassMin( double m ) {
  updated = false;
  massMin = m;
//  massSel->setMassMin( m );
//    dzSel->setMassMin( m );
  return;
}


void BPHPhiToKKBuilder::setMassMax( double m ) {
  updated = false;
  massMax = m;
//  massSel->setMassMax( m );
//    dzSel->setMassMax( m );
  return;
}


void BPHPhiToKKBuilder::setProbMin( double p ) {
  updated = false;
  probMin = p;
//  chi2Sel->setProbMin( p );
  return;
}


// mass constraint: foreseen but actually not implemented
void BPHPhiToKKBuilder::setConstr( double mass, double sigma ) {
  updated = false;
  cMass  = mass;
  cSigma = sigma;
  return;
}

/// get current cuts
double BPHPhiToKKBuilder::getPtMin() const {
  return ptMin;
//  return ptSel->getPtMin();
}


double BPHPhiToKKBuilder::getEtaMax() const {
  return etaMax;
//  return etaSel->getEtaMax();
}


double BPHPhiToKKBuilder::getMassMin() const {
  return massMin;
//  return massSel->getMassMin();
}


double BPHPhiToKKBuilder::getMassMax() const {
  return massMax;
//  return massSel->getMassMax();
}


double BPHPhiToKKBuilder::getProbMin() const {
  return probMin;
//  return chi2Sel->getProbMin();
}


// mass constraint: foreseen but actually not implemented
double BPHPhiToKKBuilder::getConstrMass() const {
  return cMass;
}


double BPHPhiToKKBuilder::getConstrSigma() const {
  return cSigma;
}

