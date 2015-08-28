#include "BPHAnalysis/RecoDecay/test/stubs/TestBPHRecoDecay.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "BPHAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "BPHAnalysis/RecoDecay/interface/BPHRecoSelect.h"
#include "BPHAnalysis/RecoDecay/interface/BPHRecoCandidate.h"
#include "BPHAnalysis/RecoDecay/interface/BPHPlusMinusCandidate.h"
#include "BPHAnalysis/RecoDecay/interface/BPHMomentumSelect.h"
#include "BPHAnalysis/RecoDecay/interface/BPHVertexSelect.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include <iostream>
#include <string>

using namespace std;

TestBPHRecoDecay::TestBPHRecoDecay( const edm::ParameterSet& ps ) {
}

TestBPHRecoDecay::~TestBPHRecoDecay() {
}

void TestBPHRecoDecay::beginJob() {
  cout << "TestBPHRecoDecay::beginJob" << endl;
  return;
}

void TestBPHRecoDecay::analyze( const edm::Event& ev,
                                const edm::EventSetup& es ) {

  cout << "--------- event "
       << ev.id().run() << " / "
       << ev.id().event() << " ---------" << endl;

//  if ( ev.id().event() != 3670523273 ) return;

  // get object collections

  edm::Handle<pat::MuonCollection> muons;
  string muonLabel = "calibratedPatMuonsPFlow";
  edm::InputTag muonTag( muonLabel );
  ev.getByLabel( muonTag, muons );
  if ( muons.isValid() ) cout << muons->size() << " muons found" << endl;
  else                   cout << "no muons" << endl;

  edm::Handle< vector<reco::Track> > tracks;
  string trackLabel = "generalTracks";
  edm::InputTag trackTag( trackLabel );
  ev.getByLabel( trackTag, tracks );
  if ( tracks.isValid() ) cout << tracks->size() << " tracks found" << endl;
  else                    cout << "no tracks" << endl;

  edm::Handle< vector<reco::PFCandidate> > pfcands;
  string pfcandLabel = "selectedPatJetsPFlow:pfCandidates:PAT";
  edm::InputTag pftag( pfcandLabel );
  ev.getByLabel( pftag, pfcands );
  if ( pfcands.isValid() ) cout << pfcands->size() << " pfcands found" << endl;
  else                     cout << "no pfcands" << endl;

  //
  // starting objects selection
  //

  // muon selection by charge
  class MuonChargeSelect: public BPHRecoSelect {
   public:
    MuonChargeSelect( int c ): charge ( c ) {}
    ~MuonChargeSelect() {}
    virtual bool accept( const reco::Candidate& cand ) const {
      const pat::Muon* p = reinterpret_cast<const pat::Muon*>( &cand );
      if ( p == 0 ) return false;
      return ( ( charge * cand.charge() ) > 0 );
    }
   private:
    int charge;
  };

  // muon selection by Pt
  class MuonPtSelect: public BPHRecoSelect {
   public:
    MuonPtSelect( float pt ): ptCut( pt ) {}
    ~MuonPtSelect() {}
    virtual bool accept( const reco::Candidate& cand ) const {
      const pat::Muon* p = reinterpret_cast<const pat::Muon*>( &cand );
      if ( p == 0 ) return false;
      return ( p->p4().pt() > ptCut );
    }
   private:
    float ptCut;
  };

  // muon selection by eta
  class MuonEtaSelect: public BPHRecoSelect {
   public:
    MuonEtaSelect( float eta ): etaCut( eta ) {}
    ~MuonEtaSelect() {}
    virtual bool accept( const reco::Candidate& cand ) const {
      const pat::Muon* p = reinterpret_cast<const pat::Muon*>( &cand );
      if ( p == 0 ) return false;
      return ( fabs( p->p4().eta() ) < etaCut );
    }
   private:
    float etaCut;
  };

  // kaon selection by charge
  class KaonChargeSelect: public BPHRecoSelect {
   public:
    KaonChargeSelect( int c ): charge ( c ) {}
    ~KaonChargeSelect() {}
    virtual bool accept( const reco::Candidate& cand ) const {
      return ( ( charge * cand.charge() ) > 0 );
    }
   private:
    int charge;
  };

  class KaonNeutralVeto: public BPHRecoSelect {
   public:
    KaonNeutralVeto() {}
    ~KaonNeutralVeto() {}
    virtual bool accept( const reco::Candidate& cand ) const {
      return lround( fabs( cand.charge() ) );
    }
  };

  // kaon selection by Pt
  class KaonPtSelect: public BPHRecoSelect {
   public:
    KaonPtSelect( float pt ): ptCut( pt ) {}
    ~KaonPtSelect() {}
    virtual bool accept( const reco::Candidate& cand ) const {
      return ( cand.p4().pt() > ptCut );
    }
   private:
    float ptCut;
  };

  // kaon selection by eta
  class KaonEtaSelect: public BPHRecoSelect {
   public:
    KaonEtaSelect( float eta ): etaCut( eta ) {}
    ~KaonEtaSelect() {}
    virtual bool accept( const reco::Candidate& cand ) const {
      return ( fabs( cand.p4().eta() ) < etaCut );
    }
   private:
    float etaCut;
  };


  //
  // reconstructed object selection
  //

  // selection by mass
  class MassSelect: public BPHMomentumSelect {
   public:
    MassSelect( double minMass, double maxMass ):
      mMin( minMass ),
      mMax( maxMass ) {}
    virtual bool accept( const BPHDecayMomentum& cand ) const {
      double mass = cand.composite().mass();
      return ( ( mass > mMin ) && ( mass < mMax ) );
    }
   private:
    double mMin;
    double mMax;
  };

  // selection by chi^2
  class Chi2Select: public BPHVertexSelect {
   public:
    Chi2Select( double minProb ):
      mProb( minProb ) {}
    virtual bool accept( const BPHDecayVertex& cand ) const {
      const reco::Vertex& v = cand.vertex();
      if ( v.isFake() ) return false;
      if ( !v.isValid() ) return false;
      return ( TMath::Prob( v.chi2(), lround( v.ndof() ) ) > mProb );
    }
   private:
    double mProb;
  };

  // build and dump JPsi

  cout << "build and dump JPsi" << endl;
  MuonPtSelect     muPt ( 4.0 );
  MuonEtaSelect    muEta( 2.1 );
  string muPos = "muPos";
  string muNeg = "muNeg";
  BPHRecoBuilder bJPsi( es );
  bJPsi.add( muPos, muons );
  bJPsi.add( muNeg, muons );
  bJPsi.filter( muPos, muPt  );
  bJPsi.filter( muNeg, muPt  );
  bJPsi.filter( muPos, muEta );
  bJPsi.filter( muNeg, muEta );

  MassSelect massJPsi( 2.5, 3.7 );
  Chi2Select chi2Valid( 0.0 );
  bJPsi.filter( massJPsi );
  bJPsi.filter( chi2Valid );
  vector<const BPHPlusMinusCandidate*> lJPsi =
               BPHPlusMinusCandidate::build( bJPsi, muPos, muNeg,
                                             3.096916, 0.00004 );
//  //  BPHPlusMinusCandidate::build function has embedded charge selection
//  //  alternatively simple BPHRecoCandidate::build function can be used
//  //  as in the following
//  MuonChargeSelect mqPos( +1 );
//  MuonChargeSelect mqNeg( +1 );
//  bJPsi.filter( nPos, mqPos );
//  bJPsi.filter( nNeg, mqNeg );
//  vector<const BPHRecoCandidate*> lJPsi =
//               BPHRecoCandidate::build( bJPsi, 
//                                        3.096916, 0.00004 );
  int iJPsi;
  int nJPsi = lJPsi.size();
  cout << nJPsi << " JPsi cand found" << endl;
  for ( iJPsi = 0; iJPsi < nJPsi; ++iJPsi ) dumpRecoCand( "JPsi",
                                                          lJPsi[iJPsi] );

  // build and dump Phi

  cout << "build and dump Phi" << endl;
  BPHRecoBuilder bPhi( es );
  KaonChargeSelect tkPos( +1 );
  KaonChargeSelect tkNeg( -1 );
  KaonPtSelect tkPt( 0.7 );
  string kPos = "kPos";
  string kNeg = "kNeg";
  bPhi.add( kPos, pfcands, 0.493677 );
  bPhi.add( kNeg, pfcands, 0.493677 );
  bPhi.filter( kPos, tkPos );
  bPhi.filter( kNeg, tkNeg );
  bPhi.filter( kPos, tkPt  );
  bPhi.filter( kNeg, tkPt  );

  MassSelect massPhi( 1.00, 1.04 );
  bPhi.filter( massPhi );
  vector<const BPHRecoCandidate*> lPhi =
               BPHRecoCandidate::build( bPhi );
//  //  BPHRecoCandidate::build function requires explicit charge selection
//  //  alternatively BPHPlusMinusCandidate::build function can be used
//  //  as in the following
//  //  (filter functions with tkPos and tkNeg can be dropped)
//  vector<const BPHPlusMinusCandidate*> lPhi =
//               BPHPlusMinusCandidate::build( bPhi, kPos, kNeg );
  int iPhi;
  int nPhi = lPhi.size();
  cout << nPhi << " Phi cand found" << endl;
  for ( iPhi = 0; iPhi < nPhi; ++iPhi ) dumpRecoCand( "Phi",
                                                      lPhi[iPhi] );

  // build and dump Bs

  if ( nJPsi && nPhi ) {
  cout << "build and dump Bs" << endl;
  BPHRecoBuilder bBs( es );
  bBs.setMinPDiffererence( 1.0e-5 );
  bBs.add( "JPsi", lJPsi );
  bBs.add(  "Phi",  lPhi );
  MassSelect mJPsi( 2.946916, 3.246916 );
  MassSelect  mPhi( 1.009461, 1.029461 );
  bBs.filter( "JPsi", mJPsi );
  bBs.filter(  "Phi",  mPhi );
  Chi2Select chi2Bs( 0.02 );
  bBs.filter( chi2Bs );
  std::vector<const BPHRecoCandidate*> lBs = BPHRecoCandidate::build( bBs );
  int iBs;
  int nBs = lBs.size();
  cout << nBs << " Bs cand found" << endl;
  for ( iBs = 0; iBs < nBs; ++iBs ) dumpRecoCand( "Bs",
                                                  lBs[iBs] );
  }

  // build and dump B+

  if ( nJPsi && pfcands->size() ) {
  cout << "build and dump B+" << endl;
  BPHRecoBuilder bBp( es );
  bBp.setMinPDiffererence( 1.0e-5 );
  bBp.add( "JPsi", lJPsi );
  bBp.add( "Kaon", pfcands, 0.493677 );
  MassSelect mJPsi( 2.946916, 3.246916 );
  KaonNeutralVeto knv;
  bBp.filter( "JPsi", mJPsi );
  bBp.filter( "Kaon", tkPt );
  bBp.filter( "Kaon", knv );
  Chi2Select chi2Bp( 0.02 );
  bBp.filter( chi2Bp );
  std::vector<const BPHRecoCandidate*> lBp = BPHRecoCandidate::build( bBp );
  int iBp;
  int nBp = lBp.size();
  cout << nBp << " B+ cand found" << endl;
  for ( iBp = 0; iBp < nBp; ++iBp ) dumpRecoCand( "B+",
                                                  lBp[iBp] );
  }

  BPHRecoCandidate::clear();

  return;

}


void TestBPHRecoDecay::endJob() {
  cout << "TestBPHRecoDecay::endJob" << endl;
  return;
}


void TestBPHRecoDecay::dumpRecoCand( const string& name,
                                     const BPHRecoCandidate* cand ) {

  static string cType = " cowboy";
  static string sType = " sailor";
  static string dType = "";
  string* type;
  const BPHPlusMinusCandidate* pmCand =
        dynamic_cast<const BPHPlusMinusCandidate*>( cand );
  if ( pmCand != 0 ) {
    if ( pmCand->isCowboy() ) type = &cType;
    else                      type = &sType;
  }
  else                        type = &dType;

  bool constrMass = ( cand->constrMass() > 0.0 );

  cout << "****** " << name << "   cand mass: "
       << cand->composite().mass() << " momentum "
       << cand->composite().px() << " "
       << cand->composite().py() << " "
       << cand->composite().pz() << *type << endl;

//  if ( constrMass )
//  cout << "       " << name << " constr mass: "
  cout << "****** " << name << " constr mass: "
       << cand->p4().mass() << " momentum "
       << cand->p4().px() << " "
       << cand->p4().py() << " "
       << cand->p4().pz() << endl;

  const reco::Vertex& vx = cand->vertex();
  const reco::Vertex::Point& vp = vx.position();
  double chi2 = vx.chi2();
  int    ndof = lround( vx.ndof() );
  double prob = TMath::Prob( chi2, ndof );
  string tdca = "";
  if ( pmCand != 0 ) {
    stringstream sstr;
    sstr << " - " << pmCand->cAppInRPhi().distance();
    tdca = sstr.str();
  }
  cout << "****** " << name << " vertex: "
       << vx.isFake() << " " << vx.isValid() << " - "
       << chi2 << " " << ndof << " " << prob << " - "
       << vp.X() << " " << vp.Y() << " " << vp.Z() << tdca << endl;

  const vector<string>& dl = cand->daugNames();
  int i;
  int n = dl.size();
  for ( i = 0; i < n; ++i ) {
    const string& name = dl[i];
    const reco::Candidate* dp = cand->getDaug( name );
    const reco::TransientTrack& tt = *cand->getTransientTrack( dp );
    GlobalPoint gp( vp.X(), vp.Y(), vp.Z() ); 
    TrajectoryStateClosestToPoint tscp =
                                  tt.trajectoryStateClosestToPoint( gp );
    GlobalVector dm = tscp.momentum();
//    TrajectoryStateOnSurface tsos = tt->stateOnSurface( gp );
//    GlobalVector gv = tsos.globalMomentum();
    cout << "daughter " << i
         << " " << name
         << " momentum: "
         << dp->px() << " "
         << dp->py() << " "
         << dp->pz() << " - at vertex: "
         << dm.x() << " "
         << dm.y() << " "
         << dm.z() << endl;
  }
  const vector<string>& dc = cand->compNames();
  int j;
  int m = dc.size();
  for ( j = 0; j < m; ++j ) {
    const string& name = dc[j];
    const BPHRecoCandidate* dp = cand->getComp( name );
    cout << "composite daughter " << j
         << " " << name
         << " momentum: "
         << dp->composite().px() << " "
         << dp->composite().py() << " "
         << dp->composite().pz() << endl;
  }

  if ( constrMass ) {
  const RefCountedKinematicTree& kt = cand->kinematicTree();
  const RefCountedKinematicVertex& kd = kt->currentDecayVertex();
  GlobalPoint gp = kd->position(); 
  cout << "   kin fit vertex: "
       << gp.x() << " "
       << gp.y() << " "
       << gp.z() << endl;
  const KinematicState& ks = kt->currentParticle()->currentState();
  GlobalVector gv = ks.globalMomentum();
  cout << "   kin fit momentum: "
       << ks.mass() << " - "
       << gv.x() << " "
       << gv.y() << " "
       << gv.z() << " - deltaM: "
       << ks.mass() - cand->constrMass() << endl;
  vector<RefCountedKinematicParticle> dk = kt->finalStateParticles();
  int j;
  int m = dk.size();
  for ( j = 0; j < m; ++j ) {
    const reco::TransientTrack& tt = dk[j]->refittedTransientTrack();
    TrajectoryStateClosestToPoint tscp =
                                  tt.trajectoryStateClosestToPoint( gp );
    GlobalVector dm = tscp.momentum();
//    TrajectoryStateOnSurface tsos = tt.stateOnSurface( gp );
//    GlobalVector dm = tsos.globalMomentum();
    cout << "daughter " << j << " refitted: "
         << dm.x() << " "
         << dm.y() << " "
         << dm.z() << endl;
  }
  }

  return;

}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( TestBPHRecoDecay );
