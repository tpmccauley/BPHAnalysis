#include "BPHAnalysis/SpecificDecay/test/stubs/CheckBPHWriteDecay.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "BPHAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "BPHAnalysis/RecoDecay/interface/BPHRecoSelect.h"
#include "BPHAnalysis/RecoDecay/interface/BPHRecoCandidate.h"
#include "BPHAnalysis/RecoDecay/interface/BPHPlusMinusCandidate.h"
#include "BPHAnalysis/RecoDecay/interface/BPHMomentumSelect.h"
#include "BPHAnalysis/RecoDecay/interface/BPHVertexSelect.h"
#include "BPHAnalysis/RecoDecay/interface/BPHTrackReference.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicState.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include <TH1.h>
#include <TFile.h>

#include <set>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

//#define SET_LABEL(NAME,PSET) ( NAME = getParameter( PSET, #NAME ) )
//// SET_LABEL(xyz,ps);
//// is equivalent to
//// xyz = getParameter( ps, "xyx" )

CheckBPHWriteDecay::CheckBPHWriteDecay( const edm::ParameterSet& ps ) {

  if ( ps.exists( "runNumber" ) )
                   runNumber = ps.getParameter<unsigned int>( "runNumber" );
  else             runNumber = 0;
  if ( ps.exists( "evtNumber" ) )
                   evtNumber = ps.getParameter<unsigned int>( "evtNumber" );
  else             evtNumber = 0;
/*
  SET_LABEL( oniaCandsLabel, ps );
  SET_LABEL(   sdCandsLabel, ps );
  SET_LABEL(   ssCandsLabel, ps );
  SET_LABEL(   buCandsLabel, ps );
  SET_LABEL(   bdCandsLabel, ps );
  SET_LABEL(   bsCandsLabel, ps );
  consume< vector<pat::CompositeCandidate> >( oniaCandsToken, oniaCandsLabel );
  consume< vector<pat::CompositeCandidate> >(   sdCandsToken,   sdCandsLabel );
  consume< vector<pat::CompositeCandidate> >(   ssCandsToken,   ssCandsLabel );
  consume< vector<pat::CompositeCandidate> >(   buCandsToken,   buCandsLabel );
  consume< vector<pat::CompositeCandidate> >(   bdCandsToken,   bdCandsLabel );
  consume< vector<pat::CompositeCandidate> >(   bsCandsToken,   bsCandsLabel );
*/
  candsLabel = ps.getParameter< std::vector<std::string> >( "candsLabel" );
  int i;
  int n =
  candsLabel.  size();
  candsToken.resize( n );
  for ( i = 0; i < n; ++i )
        consume< vector<pat::CompositeCandidate> >( candsToken[i],
                                                    candsLabel[i] );

  if ( ps.exists( "fileName" ) )
       osPtr = new ofstream( ps.getParameter<string>( "fileName" ) );
  else osPtr = &cout;


}


CheckBPHWriteDecay::~CheckBPHWriteDecay() {
}


void CheckBPHWriteDecay::beginJob() {
  return;
}


void CheckBPHWriteDecay::analyze( const edm::Event& ev,
                                  const edm::EventSetup& es ) {

  if ( ( runNumber != 0 ) && ( ev.id().run  () != runNumber ) ) return;
  if ( ( evtNumber != 0 ) && ( ev.id().event() != evtNumber ) ) return;
  cout << "--------- event "
       << ev.id().run() << " / "
       << ev.id().event() << " ---------" << endl;

  ostream& os = *osPtr;

  int il;
  int nl =
  candsLabel.size();
  for ( il = 0; il < nl; ++il ) {
    edm::Handle< vector<pat::CompositeCandidate> > cands;
    candsToken[il].get( ev, cands );
    int ic;
    int nc = cands->size();
    for ( ic = 0; ic < nc; ++ ic ) {
      os << "*********** " << candsLabel[il] << " " << ic << "/" << nc
         << " ***********"
         << endl;
      const pat::CompositeCandidate& cand = cands->at( ic );
      dump( os, cand );
    }
  }
/*
  //////////// quarkonia ////////////

  edm::Handle< vector<pat::CompositeCandidate> > oniaCands;
  oniaCandsToken.get( ev, oniaCands );
  int iqo;
  int nqo = oniaCands->size();
  for ( iqo = 0; iqo < nqo; ++ iqo ) {
    cout << "*********** quarkonium " << iqo << "/" << nqo << " ***********"
         << endl;
    const pat::CompositeCandidate& cand = oniaCands->at( iqo );
    dump ( cand );
  }
  os << "onia " << nqo << endl;

  //////////// Bu ////////////

  edm::Handle< vector<pat::CompositeCandidate> > buCands;
  buCandsToken.get( ev, buCands );
  int ibu;
  int nbu = buCands->size();
  for ( ibu = 0; ibu < nbu; ++ ibu ) {
    cout << "*********** Bu " << ibu << "/" << nbu << " ***********"
         << endl;
    const pat::CompositeCandidate& cand = buCands->at( ibu );
    dump ( cand );
  }
  os << "Bu " << nbu << endl;

  //////////// Bd ////////////

  edm::Handle< vector<pat::CompositeCandidate> > bdCands;
  bdCandsToken.get( ev, bdCands );
  int ibd;
  int nbd = bdCands->size();
  for ( ibd = 0; ibd < nbd; ++ ibd ) {
    cout << "*********** Bd " << ibd << "/" << nbd << " ***********"
         << endl;
    const pat::CompositeCandidate& cand = bdCands->at( ibd );
    dump ( cand );
  }
  os << "Bd " << nbd << endl;

  //////////// Sd ////////////

  edm::Handle< vector<pat::CompositeCandidate> > sdCands;
  sdCandsToken.get( ev, sdCands );
  int isd;
  int nsd = sdCands->size();
  for ( isd = 0; isd < nsd; ++ isd ) {
    cout << "*********** Sd " << isd << "/" << nsd << " ***********"
         << endl;
    const pat::CompositeCandidate& cand = sdCands->at( isd );
    dump ( cand );
  }
  os << "Sd " << nsd << endl;

  //////////// Bs ////////////

  edm::Handle< vector<pat::CompositeCandidate> > bsCands;
  bsCandsToken.get( ev, bsCands );
  int ibs;
  int nbs = bsCands->size();
  for ( ibs = 0; ibs < nbs; ++ ibs ) {
    cout << "*********** Bs " << ibs << "/" << nbs << " ***********"
         << endl;
    const pat::CompositeCandidate& cand = bsCands->at( ibs );
    dump ( cand );
  }
  os << "Bs " << nbs << endl;

  //////////// Ss ////////////

  edm::Handle< vector<pat::CompositeCandidate> > ssCands;
  ssCandsToken.get( ev, ssCands );
  int iss;
  int nss = ssCands->size();
  for ( iss = 0; iss < nss; ++ iss ) {
    cout << "*********** Ss " << iss << "/" << nss << " ***********"
         << endl;
    const pat::CompositeCandidate& cand = ssCands->at( iss );
    dump ( cand );
  }
  os << "Ss " << nss << endl;
*/
  return;

}


void CheckBPHWriteDecay::endJob() {
  return;
}


void CheckBPHWriteDecay::dump( std::ostream& os,
                               const pat::CompositeCandidate& cand ) {

  float mfit = ( cand.hasUserFloat( "fitMass" ) ?
                 cand.   userFloat( "fitMass" ) : -1 );
  os << &cand
     << " mass : " << cand.mass() << " " << mfit << " "
     << (   cand.hasUserData      ( "cowboy" ) ?
          ( cand.   userData<bool>( "cowboy" ) ? "cowboy" : "sailor" )
                                                   : "" ) << endl;
  writeMomentum( os, "cmom ", cand, false );
  writePosition( os, " xyz ", cand.momentum() );
  const reco::Vertex* vptr = 
        ( cand.hasUserData              ( "vertex" ) ? 
          cand.   userData<reco::Vertex>( "vertex" ) : 0 );
  if ( vptr != 0 ) {
    writePosition( os, "vpos : ", *vptr, false );
    os << " --- " << vptr->chi2() << " / " << vptr->ndof()
       << " ( " << ChiSquaredProbability( vptr->chi2(),
                                            vptr->ndof() ) << " ) " << endl;
  }
  const reco::Vertex* vfit = 
        ( cand.hasUserData              ( "fitVertex" ) ? 
          cand.   userData<reco::Vertex>( "fitVertex" ) : 0 );
  if ( vfit != 0 ) {
    writePosition( os, "vfit : ", *vfit, false );
    os << " --- "  << vfit->chi2() << " / " << vfit->ndof()
       << " ( " << ChiSquaredProbability( vfit->chi2(),
                                            vfit->ndof() ) << " ) " << endl;
  }
  if ( cand.hasUserData( "fitMomentum" ) )
       writePosition( os, "fmom : ",
      *cand.   userData< Vector3DBase<float,GlobalTag> >( "fitMomentum" ) );

  if ( cand.hasUserData( "primaryVertex" ) ) {
    const vertex_ref* pvr = cand.userData<vertex_ref>( "primaryVertex" );
    if ( pvr->isNonnull() ) {
      const reco::Vertex* pvtx = pvr->get();
      if ( pvtx != 0 ) writePosition( os, "ppos ", *pvtx );
    }
  }
  int i;
  int n = cand.numberOfDaughters();
  for ( i = 0; i < n; ++i ) {
    const reco::Candidate* dptr = cand.daughter( i );
    os << "daug " << i << "/" << n << " : " << dptr;
    writeMomentum( os, " == ", *dptr, false );
    os << " " << dptr->mass() << endl;
    const pat::Muon* mptr = dynamic_cast<const pat::Muon*>( dptr );
    os << "muon " << i << "/" << n << " : " << mptr << endl;
    const reco::Track* tptr = BPHTrackReference::getTrack( *dptr, "cfhpmnigs" );
    os << "trk  " << i << "/" << n << " : " << tptr;
    if ( tptr != 0 ) writeMomentum( os, " == ", *tptr );
    else             os << "no track" << endl;
  }
  const vector<string>& names = cand.userDataNames();
  int j;
  int m = names.size();
  for ( j = 0; j < m; ++j ) {
    const string& dname = names[j];
    if ( dname.substr( 0, 5 ) != "refTo" ) continue;
    const compcc_ref* ref = cand.userData<compcc_ref>( dname );
    os << dname << " : " << ref->get() << endl;
  }

  return;

}


//string CheckBPHWriteDecay::getParameter( const edm::ParameterSet& ps,
//                                         const string& name ) {
//  if ( ps.exists( name ) ) return ps.getParameter<string>( name );
//  return "";
//}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( CheckBPHWriteDecay );
