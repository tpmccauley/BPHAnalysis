/*
 *  See header file for a description of this class.
 *
 *  $Date: 2015-07-03 16:33:59 $
 *  $Revision: 1.1 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//-----------------------
// This Class' Header --
//-----------------------
#include "BPHAnalysis/RecoDecay/interface/BPHPlusMinusCandidate.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BPHAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "BPHAnalysis/RecoDecay/interface/BPHRecoSelect.h"
#include "BPHAnalysis/RecoDecay/interface/BPHMomentumSelect.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"

//---------------
// C++ Headers --
//---------------
#include <iostream>
#include <math.h>

using namespace std;

//-------------------
// Initializations --
//-------------------


//----------------
// Constructors --
//----------------
BPHPlusMinusCandidate::BPHPlusMinusCandidate( const edm::EventSetup* es ):
  BPHDecayVertex( es ),
  BPHPlusMinusVertex( es ),
  BPHRecoCandidate( es ) {
}


BPHPlusMinusCandidate::BPHPlusMinusCandidate( const edm::EventSetup* es,
                       const BPHRecoBuilder::ComponentSet& compList ):
  BPHDecayMomentum( compList.daugMap, compList.compMap ),
  BPHDecayVertex( this, es ),
  BPHKinematicFit( this ),
  BPHPlusMinusVertex( es ),
  BPHRecoCandidate( es, compList ) {
}

//--------------
// Destructor --
//--------------
BPHPlusMinusCandidate::~BPHPlusMinusCandidate() {
}

//--------------
// Operations --
//--------------
void BPHPlusMinusCandidate::add( const std::string& name,
                                 const reco::Candidate* daug, double mass ) {
  const vector<const reco::Candidate*>& dL = daughters();
  switch ( dL.size() ) {
  case 2:
    cout << "BPHPlusMinusCandidate complete, add rejected" << endl;
    return;
  case 1:
    if ( ( daug->charge() * dL.front()->charge() ) > 0 ) {
      cout << "BPHPlusMinusCandidate already containing same sign "
           << "particle, add rejected" << endl;
      return;
    }
  case 0:
    BPHDecayMomentum::add( name, daug );
  }
  return;
}


std::vector<BPHPlusMinusConstCandPtr> BPHPlusMinusCandidate::build(
                                      const BPHRecoBuilder& builder,
                                      const std::string nPos,
                                      const std::string nNeg,
                                      double mass, double msig ) {
  vector<BPHPlusMinusConstCandPtr> cList;
  class ChargeSelect: public BPHRecoSelect {
   public:
    ChargeSelect( int c ): charge ( c ) {}
    ~ChargeSelect() {}
    virtual bool accept( const reco::Candidate& cand ) const {
      return ( ( charge * cand.charge() ) > 0 );
    }
   private:
    int charge;
  };
  ChargeSelect tkPos( +1 );
  ChargeSelect tkNeg( -1 );
  builder.filter( nPos, tkPos );
  builder.filter( nNeg, tkNeg );
  fill( cList, builder, mass, msig );
  return cList;
}


const pat::CompositeCandidate& BPHPlusMinusCandidate::composite() const {
  static pat::CompositeCandidate compCand;
  static string msg =
         "BPHPlusMinusCandidate incomplete, no composite available";
  if ( !chkSize( msg ) ) {
    compCand.clearDaughters();
    return compCand;
  }
  return BPHDecayMomentum::composite();
}


bool BPHPlusMinusCandidate::isCowboy() const {
  static string msg =
         "BPHPlusMinusCandidate incomplete, no cowboy/sailor classification";
  return ( chkSize( msg ) && phiDiff() );
}


bool BPHPlusMinusCandidate::isSailor() const {
  static string msg =
         "BPHPlusMinusCandidate incomplete, no cowboy/sailor classification";
  return ( chkSize( msg ) && !phiDiff() );
}


bool BPHPlusMinusCandidate::phiDiff() const {
  const vector<const reco::Candidate*>& dL = daughters();
  int idPos = ( dL.front()->charge() > 0 ? 0 : 1 );
  double dphi = dL[idPos]->phi() - dL[1 - idPos]->phi();
  double TWO_PI = 2 * M_PI;
  while ( dphi >  M_PI ) dphi -= TWO_PI;
  while ( dphi < -M_PI ) dphi += TWO_PI;
  return dphi > 0;
}

