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
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
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
  evSetup( es ),
  updatedTracks( false ),
  updatedVertex( false ) {
}

//--------------
// Destructor --
//--------------
BPHDecayVertex::~BPHDecayVertex() {
}

//--------------
// Operations --
//--------------
const reco::Vertex& BPHDecayVertex::vertex() const {
  fitVertex();
  return fittedVertex;
}


const vector<reco::TransientTrack>& BPHDecayVertex::transientTracks() const {
  tTracks();
  return trTracks;
}


reco::TransientTrack* BPHDecayVertex::getTransientTrack(
                                      const reco::Candidate* cand ) const {
  tTracks();
  map<const reco::Candidate*,
            reco::TransientTrack*>::const_iterator iter = ttMap.find( cand );
  return ( iter != ttMap.end() ? iter->second : 0 );
}


void BPHDecayVertex::setNotUpdated() const {
  BPHDecayMomentum::setNotUpdated();
  updatedTracks = updatedVertex = false;
  return;
}


void BPHDecayVertex::tTracks() const {
  if ( updatedTracks ) return;
  trTracks.clear();
  ttMap.clear();
  edm::ESHandle<TransientTrackBuilder> ttB;
  evSetup->get<TransientTrackRecord>().get( "TransientTrackBuilder", ttB );
  const vector<const reco::Candidate*>& dL = daughFull();
  int n = dL.size();
  trTracks.reserve( n );
  while ( n-- ) {
    const reco::Candidate* rp = dL[n];
    if ( !rp->charge() ) {
      ttMap[rp] = 0;
      continue;
    }
    const reco::Track* tp;
    const reco::PFCandidate* pp =
          dynamic_cast<const reco::PFCandidate*>( rp );
    if ( pp != 0 ) tp = pp->trackRef().get();
    else           tp = rp->get<reco::TrackRef>().get();
    if ( tp == 0 ) {
      cout << "no track" << endl;
      continue;
    }
    trTracks.push_back( ttB->build( tp ) );
    reco::TransientTrack* ttp = &trTracks.back();
    ttMap[rp] = ttp;
  }
  updatedTracks = true;
  return;
}


void BPHDecayVertex::fitVertex() const {
  if ( updatedVertex ) return;
  tTracks();
  KalmanVertexFitter kvf( true );
  TransientVertex tv = kvf.vertex( trTracks );
  fittedVertex = tv;
  updatedVertex = true;
  return;
}

