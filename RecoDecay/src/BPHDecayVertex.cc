/*
 *  See header file for a description of this class.
 *
 *  $Date: 2015-07-03 13:49:53 $
 *  $Revision: 1.1 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//-----------------------
// This Class' Header --
//-----------------------
#include "BPHAnalysis/RecoDecay/interface/BPHDecayVertex.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BPHAnalysis/RecoDecay/interface/BPHRecoCandidate.h"
#include "BPHAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "BPHAnalysis/RecoDecay/interface/BPHTrackReference.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
//#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

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
BPHDecayVertex::BPHDecayVertex( const edm::EventSetup* es ):
 evSetup( es ) {
  setNotUpdated();
}

//--------------
// Destructor --
//--------------
BPHDecayVertex::~BPHDecayVertex() {
}

//--------------
// Operations --
//--------------
bool BPHDecayVertex::isValid() const {
  if ( oldVertex ) fitVertex();
  return validVertex;
}


const reco::Vertex& BPHDecayVertex::vertex() const {
  if ( oldVertex ) fitVertex();
  return fittedVertex;
}


const vector<reco::TransientTrack>& BPHDecayVertex::transientTracks() const {
  if ( oldTracks ) tTracks();
  return trTracks;
}


reco::TransientTrack* BPHDecayVertex::getTransientTrack(
                                      const reco::Candidate* cand ) const {
  if ( oldTracks ) tTracks();
  map<const reco::Candidate*,
            reco::TransientTrack*>::const_iterator iter = ttMap.find( cand );
  return ( iter != ttMap.end() ? iter->second : 0 );
}


void BPHDecayVertex::setNotUpdated() const {
  BPHDecayMomentum::setNotUpdated();
  oldTracks = oldVertex = true;
  validVertex = false;
  return;
}


void BPHDecayVertex::tTracks() const {
  oldTracks = false;
  trTracks.clear();
  ttMap.clear();
  edm::ESHandle<TransientTrackBuilder> ttB;
  evSetup->get<TransientTrackRecord>().get( "TransientTrackBuilder", ttB );
  const vector<const reco::Candidate*>& dL = daughFull();
  int n = dL.size();
  trTracks.reserve( n );
  validVertex = true;
  while ( n-- ) {
    const reco::Candidate* rp = dL[n];
    ttMap[rp] = 0;
    if ( !rp->charge() ) continue;
    const reco::Track* tp;
    tp = BPHTrackReference::getTrack( *rp, "cfp" );
    if ( tp == 0 ) {
      cout << "no track for reco::(PF)Candidate" << endl;
      validVertex = false;
      continue;
    }
    trTracks.push_back( ttB->build( tp ) );
    reco::TransientTrack* ttp = &trTracks.back();
    ttMap[rp] = ttp;
  }
  return;
}


void BPHDecayVertex::fitVertex() const {
  oldVertex = false;
  if ( oldTracks ) tTracks();
  if ( trTracks.size() < 2 ) return;
  KalmanVertexFitter kvf( true );
  TransientVertex tv = kvf.vertex( trTracks );
  fittedVertex = tv;
  return;
}

