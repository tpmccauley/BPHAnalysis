
#include "BPHAnalysis/SpecificDecay/plugins/BPHHistoSpecificDecay.h"

#include "BPHAnalysis/SpecificDecay/interface/BPHParticleMasses.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"

#include <TH1.h>
#include <TFile.h>

#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

#define SET_LABEL(NAME,PSET) ( NAME = PSET.getParameter<string>( #NAME ) )
// SET_LABEL(xyz,ps);
// is equivalent to
// xyz = ps.getParameter<string>( "xyx" );

#define CAT3(A,B,C) A ## B ## C
#define STRING_NX(A) #A
#define STRING(A) STRING_NX(A)
#define CHK_TRIG(RESULTS,NAMES,INDEX,PATH) if ( NAMES[INDEX].find( STRING( CAT3( HLT_, PATH, _v ) ) ) == 0 ) { flag_ ## PATH = RESULTS->accept( INDEX ); continue; }
// CHK_TRIG( trigResults, names, i, xyz );
// is equivalent to
// if ( names[i].find( "HLT_xyz_v" ) == 0 ) { flag_xyz = trigResults->accept( i ); break; }

class BPHUserData {
 public:
  template<class T>
  static const T* get( const pat::CompositeCandidate& cand,
                       const string& name ) {
    if ( cand.hasUserData( name ) ) return cand.userData<T>( name );
    return 0;
  }
  template<class T>
  static const T* getByRef( const pat::CompositeCandidate& cand,
                            const string& name ) {
    if ( cand.hasUserData( name ) ) {
      typedef edm::Ref< std::vector<T> > objRef;
      const objRef* ref = cand.userData<objRef>( name );
      if ( ref ==      0 ) return 0;
      if ( ref->isNull() ) return 0;
      return ref->get();
    }
    return 0;
  }
};


class BPHDaughters {
 public:
  static vector<const reco::Candidate*> get(
         const pat::CompositeCandidate& cand, float massMin, float massMax ) {
    int i;
    int n = cand.numberOfDaughters();
    vector<const reco::Candidate*> v;
    v.reserve( n );
    for ( i = 0; i < n; ++i ) {
      const reco::Candidate* dptr = cand.daughter( i );
      float mass = dptr->mass();
      if ( ( mass > massMin ) && ( mass < massMax ) ) v.push_back( dptr );
    }
    return v;
  }
};


class BPHSoftMuonSelect {
 public:
  BPHSoftMuonSelect( int   cutTrackerLayers = 5,
                     int   cutPixelLayers   = 0,
                     float maxDxy           = 0.3,
                     float maxDz            = 20.0,
                     bool  goodMuon         = true ,
                     bool  highPurity       = true ):
   cutTL ( cutTrackerLayers ),
   cutPL ( cutPixelLayers   ),
   maxXY ( maxDxy           ),
   maxZ  ( maxDz            ),
   gM    ( goodMuon         ),
   hP    ( highPurity       ) {}
  bool accept( const reco::Candidate& cand,
               const reco::Vertex* pv ) const {
    const pat::Muon* p = dynamic_cast<const pat::Muon*>( &cand );
    if ( p == 0 ) return false;
    if ( gM &&
        !muon::isGoodMuon( *p, muon::TMOneStationTight ) ) return false;
    if ( p->innerTrack()->hitPattern().trackerLayersWithMeasurement() <= cutTL )
         return false;
    if ( p->innerTrack()->hitPattern().  pixelLayersWithMeasurement() <= cutPL )
         return false;
    if ( hP &&
        !p->innerTrack()->quality( reco::TrackBase::highPurity ) )
         return false;
    if ( pv == 0 ) return true;
    const reco::Vertex::Point& pos = pv->position();
    if ( fabs( p->innerTrack()->dxy( pos ) ) >= maxXY )
         return false;
    if ( fabs( p->innerTrack()->dz ( pos ) ) >= maxZ )
         return false;
    return true;
  }
 private:
  const reco::Vertex* pv;
  int cutTL;
  int cutPL;
  float maxXY;
  float maxZ;
  bool gM;
  bool hP;
};


class BPHDaughterSelect: public BPHHistoSpecificDecay::CandidateSelect {
 public:
  BPHDaughterSelect( float  ptMinLoose,
                     float  ptMinTight,
                     float etaMaxLoose,
                     float etaMaxTight,
                     const BPHSoftMuonSelect*
                           softMuonselector = 0 ): pLMin(  ptMinLoose ),
                                                   pTMin(  ptMinTight ),
                                                   eLMax( etaMaxLoose ),
                                                   eTMax( etaMaxTight ),
                                                   sms( softMuonselector ) {
  }
  bool accept( const pat::CompositeCandidate& cand,
               const reco::Vertex* pv = 0 ) const {
/*
    const reco::Candidate* dptr0 = cand.daughter( 0 );
    const reco::Candidate* dptr1 = cand.daughter( 1 );
    if ( dptr0 == 0 ) return false;
    if ( dptr1 == 0 ) return false;
    float pt0 = dptr0->pt();
    float pt1 = dptr1->pt();
    if ( ( pt0 < pLMin ) || ( pt1 < pLMin ) ) return false;
    if ( ( pt0 < pTMin ) && ( pt1 < pTMin ) ) return false;
    float eta0 = fabs( dptr0->eta() );
    float eta1 = fabs( dptr1->eta() );
    if (   ( eLMax > 0 ) &&
         ( ( eta0 > eLMax ) || ( eta1 > eLMax ) ) ) return false;
    if (   ( eTMax > 0 ) && 
         ( ( eta0 > eTMax ) && ( eta1 > eTMax ) ) ) return false;
    if ( sms != 0 ) {
      const reco::Vertex* pvtx = BPHUserData::getByRef
           <reco::Vertex>( cand, "primaryVertex" );
      if ( pvtx == 0 ) return false;
      if ( !sms->accept( *dptr0, pvtx ) ) return false;
      if ( !sms->accept( *dptr1, pvtx ) ) return false;
    }
    return true;
*/
    return accept( cand, pLMin, pTMin, eLMax, eTMax, pv, sms );
  }
  static
  bool accept( const pat::CompositeCandidate& cand,
               float  ptMinLoose,
               float  ptMinTight,
               float etaMaxLoose,
               float etaMaxTight,
               const reco::Vertex* pv = 0,
               const BPHSoftMuonSelect* softMuonselector = 0 ) {
    const reco::Candidate* dptr0 = cand.daughter( 0 );
    const reco::Candidate* dptr1 = cand.daughter( 1 );
    if ( dptr0 == 0 ) return false;
    if ( dptr1 == 0 ) return false;
    float pt0 = dptr0->pt();
    float pt1 = dptr1->pt();
    if ( ( pt0 < ptMinLoose ) || ( pt1 < ptMinLoose ) ) return false;
    if ( ( pt0 < ptMinTight ) && ( pt1 < ptMinTight ) ) return false;
    float eta0 = fabs( dptr0->eta() );
    float eta1 = fabs( dptr1->eta() );
    if (   ( etaMaxLoose > 0 ) &&
         ( ( eta0 > etaMaxLoose ) || ( eta1 > etaMaxLoose ) ) ) return false;
    if (   ( etaMaxTight > 0 ) && 
         ( ( eta0 > etaMaxTight ) && ( eta1 > etaMaxTight ) ) ) return false;
    if ( softMuonselector != 0 ) {
      const reco::Vertex* pvtx = BPHUserData::getByRef
           <reco::Vertex>( cand, "primaryVertex" );
      if ( pvtx == 0 ) return false;
      if ( !softMuonselector->accept( *dptr0, pvtx ) ) return false;
      if ( !softMuonselector->accept( *dptr1, pvtx ) ) return false;
    }    
    return true;
  }
 private:
  float pLMin;
  float pTMin;
  float eLMax;
  float eTMax;
  const BPHSoftMuonSelect* sms;
};


class BPHCompositeBasicSelect: public BPHHistoSpecificDecay::CandidateSelect {
 public:
  BPHCompositeBasicSelect( float     massMin,
                           float     massMax,
                           float       ptMin = -1.0,
                           float      etaMax = -1.0,
                           float rapidityMax = -1.0 ): mMin(     massMin ),
                                                       mMax(     massMax ),
                                                       pMin(       ptMin ),
                                                       eMax(      etaMax ),
                                                       yMax( rapidityMax ) {
  }
  bool accept( const pat::CompositeCandidate& cand,
               const reco::Vertex* pv = 0 ) const {
    if ( ( ( mMin > 0 ) && ( mMax < 0 ) ) ||
         ( ( mMin < 0 ) && ( mMax > 0 ) ) ||
         ( ( mMin > 0 ) && ( mMax > 0 ) && ( mMin < mMax ) ) ) {
      float mass = cand.mass();
      if (   mass < mMin   ) return false;
      if ( ( mMax >    0 ) &&
           ( mass > mMax ) ) return false;
    }
    if (         cand.      pt()   < pMin   ) return false;
    if ( ( eMax                    >    0 ) &&
         ( fabs( cand.     eta() ) > eMax ) ) return false;
    if ( ( yMax                    >    0 ) &&
         ( fabs( cand.rapidity() ) > yMax ) ) return false;
    return true;
  }

 private:
  float mMin;
  float mMax;
  float pMin;
  float eMax;
  float yMax;
};


class BPHFittedBasicSelect: public BPHHistoSpecificDecay::CandidateSelect {
 public:
  BPHFittedBasicSelect( float     massMin,
                        float     massMax,
                        float       ptMin = -1.0,
                        float      etaMax = -1.0,
                        float rapidityMax = -1.0 ): mMin(     massMin ),
                                                    mMax(     massMax ),
                                                    pMin(       ptMin ),
                                                    eMax(      etaMax ),
                                                    yMax( rapidityMax ) {
  }
  bool accept( const pat::CompositeCandidate& cand,
               const reco::Vertex* pv = 0 ) const {
    if ( !cand.hasUserFloat( "fitMass"     ) ) return false;
    float mass = cand.userFloat( "fitMass" );
    if ( ( ( mMin > 0 ) && ( mMax < 0 ) ) ||
         ( ( mMin < 0 ) && ( mMax > 0 ) ) ||
         ( ( mMin > 0 ) && ( mMax > 0 ) && ( mMin < mMax ) ) ) {
      if (   mass < mMin   ) return false;
      if ( ( mMax >    0 ) &&
           ( mass > mMax ) ) return false;
    }
    const Vector3DBase<float,GlobalTag>* fmom = BPHUserData::get
        < Vector3DBase<float,GlobalTag> >( cand, "fitMomentum" );
    if ( fmom == 0 ) return false;
    if ( pMin > 0 ) {
      if ( fmom->transverse() < pMin ) return false;
    }
    if ( eMax > 0 ) {
      if ( fabs( fmom->eta() ) > eMax ) return false;
    }
    if ( yMax > 0 ) {
      float x = fmom->x();
      float y = fmom->y();
      float z = fmom->z();
      float e = sqrt( ( x * x ) + ( y * y ) + ( z * z ) + ( mass * mass ) );
      float r = log( ( e + z ) / ( e - z ) ) / 2;
      if ( fabs( r ) > yMax ) return false;
    }
    return true;
  }

 private:
  float mMin;
  float mMax;
  float pMin;
  float eMax;
  float yMax;
};


class BPHVertexSelect: public BPHHistoSpecificDecay::CandidateSelect {
 public:
  BPHVertexSelect( char vType,
                   float probMin,
                   float  cosMin = -1.0,
                   float  sigMin = -1.0 ): type( vType ),
                                           pMin( probMin ),
                                           cMin(  cosMin ),
                                           sMin(  sigMin ) {
  }
  bool accept( const pat::CompositeCandidate& cand,
               const reco::Vertex* pvtx ) const {
    const reco::Vertex* svtx;
    float px;
    float py;
    float mass;
    switch ( type ) {
    case 'c':
      svtx = BPHUserData::get<reco::Vertex>( cand, "vertex" );
      px = cand.px();
      py = cand.py();
      mass = cand.mass();
      break;
    case 'f':
      svtx = BPHUserData::get<reco::Vertex>( cand, "fitVertex" );
      {
      const Vector3DBase<float,GlobalTag>* fmom = BPHUserData::get
          < Vector3DBase<float,GlobalTag> >( cand, "fitMomentum" );
      if ( fmom == 0 ) return false;
      px = fmom->x();
      py = fmom->y();
      }
      if ( !cand.hasUserFloat( "fitMass" ) ) return false;
      mass = cand.userFloat( "fitMass" );
      break;
    default:
      return false;
    }
    if ( svtx == 0 ) return false;
    if ( pvtx == 0 ) return false;
    if ( pMin > 0 ) {
      if ( ChiSquaredProbability( svtx->chi2(),
                                  svtx->ndof() ) < pMin ) return false;
    }
    if ( ( cMin > 0 ) || ( sMin > 0 ) ) {
      TVector3 disp( svtx->x() - pvtx->x(),
                     svtx->y() - pvtx->y(),
                     0 );
      TVector3 cmom( px, py, 0 );
      float cosAlpha = disp.Dot( cmom ) / ( disp.Perp() * cmom.Perp() );
      if ( cosAlpha < cMin ) return false;
      if ( sMin < 0 ) return true;
      AlgebraicVector3 vmom( px, py, 0 );
      VertexDistanceXY vdistXY;
      Measurement1D distXY = vdistXY.distance( *svtx, *pvtx );
      double ctauPV = distXY.value() * cosAlpha * mass / cmom.Perp();
      GlobalError sve = svtx->error();
      GlobalError pve = pvtx->error();
      AlgebraicSymMatrix33 vXYe = sve.matrix() + pve.matrix();
      double ctauErrPV = sqrt( ROOT::Math::Similarity( vmom, vXYe ) ) * mass /
                               cmom.Perp2();
      if ( ( ctauPV / ctauErrPV ) < sMin ) return false;
    }
    return true;
  }

 private:
  char type;
  float pMin;
  float cMin;
  float sMin;
};


class BPHCompositeVertexSelect: public BPHHistoSpecificDecay::CandidateSelect {
 public:
  BPHCompositeVertexSelect( float probMin,
                            float  cosMin = -1.0,
                            float  sigMin = -1.0 ): pMin( probMin ),
                                                    cMin(  cosMin ),
                                                    sMin(  sigMin ) {
  }
  bool accept( const pat::CompositeCandidate& cand,
               const reco::Vertex* pvtx ) const {
    const reco::Vertex* svtx = BPHUserData::get
         <reco::Vertex>( cand, "vertex" );
    if ( svtx == 0 ) return false;
    if ( pvtx == 0 ) return false;
    if ( pMin > 0 ) {
      if ( ChiSquaredProbability( svtx->chi2(),
                                  svtx->ndof() ) < pMin ) return false;
    }
    if ( ( cMin > 0 ) || ( sMin > 0 ) ) {
      TVector3 disp( svtx->x() - pvtx->x(),
                     svtx->y() - pvtx->y(),
                     0 );
      TVector3 cmom( cand.px(), cand.py(), 0 );
      float cosAlpha = disp.Dot( cmom ) / ( disp.Perp() * cmom.Perp() );
      if ( cosAlpha < cMin ) return false;
      if ( sMin < 0 ) return true;
      float mass = cand.mass();
      AlgebraicVector3 vmom( cand.px(), cand.py(), 0 );
      VertexDistanceXY vdistXY;
      Measurement1D distXY = vdistXY.distance( *svtx, *pvtx );
      double ctauPV = distXY.value() * cosAlpha * mass / cmom.Perp();
      GlobalError sve = svtx->error();
      GlobalError pve = pvtx->error();
      AlgebraicSymMatrix33 vXYe = sve.matrix() + pve.matrix();
      double ctauErrPV = sqrt( ROOT::Math::Similarity( vmom, vXYe ) ) * mass /
                               cmom.Perp2();
      if ( ( ctauPV / ctauErrPV ) < sMin ) return false;
    }
    return true;
  }

 private:
  float pMin;
  float cMin;
  float sMin;
};


class BPHFittedVertexSelect: public BPHHistoSpecificDecay::CandidateSelect {
 public:
  BPHFittedVertexSelect( float probMin,
                         float  cosMin = -1.0,
                         float  sigMin = -1.0 ): pMin( probMin ),
                                                 cMin(  cosMin ),
                                                 sMin(  sigMin ) {
  }
  bool accept( const pat::CompositeCandidate& cand,
               const reco::Vertex* pvtx ) const {
    const reco::Vertex* svtx = BPHUserData::get
         <reco::Vertex>( cand, "fitVertex" );
    if ( svtx == 0 ) return false;
    if ( pvtx == 0 ) return false;
    if ( pMin > 0 ) {
      if ( ChiSquaredProbability( svtx->chi2(),
                                  svtx->ndof() ) < pMin ) return false;
    }
    if ( ( cMin > 0 ) || ( sMin > 0 ) ) {
      TVector3 disp( svtx->x() - pvtx->x(),
                     svtx->y() - pvtx->y(),
                     0 );
      const Vector3DBase<float,GlobalTag>* fmom = BPHUserData::get
          < Vector3DBase<float,GlobalTag> >( cand, "fitMomentum" );
      if ( fmom == 0 ) return false;
      TVector3 cmom( fmom->x(), fmom->y(), 0 );
      float cosAlpha = disp.Dot( cmom ) / ( disp.Perp() * cmom.Perp() );
      if ( cosAlpha < cMin ) return false;
      if ( sMin < 0 ) return true;
      if ( !cand.hasUserFloat( "fitMass" ) ) return false;
      float mass = cand.userFloat( "fitMass" );
      AlgebraicVector3 vmom( fmom->x(), fmom->y(), 0 );
      VertexDistanceXY vdistXY;
      Measurement1D distXY = vdistXY.distance( *svtx, *pvtx );
      double ctauPV = distXY.value() * cosAlpha * mass / cmom.Perp();
      GlobalError sve = svtx->error();
      GlobalError pve = pvtx->error();
      AlgebraicSymMatrix33 vXYe = sve.matrix() + pve.matrix();
      double ctauErrPV = sqrt( ROOT::Math::Similarity( vmom, vXYe ) ) * mass /
                               cmom.Perp2();
      if ( ( ctauPV / ctauErrPV ) < sMin ) return false;
    }
    return true;
  }

 private:
  float pMin;
  float cMin;
  float sMin;
};


BPHHistoSpecificDecay::BPHHistoSpecificDecay( const edm::ParameterSet& ps ) {

  useTrig = ( SET_LABEL( trigResultsLabel, ps ) != "" );
  useOnia = ( SET_LABEL(   oniaCandsLabel, ps ) != "" );
  useSd   = ( SET_LABEL(     sdCandsLabel, ps ) != "" );
  useSs   = ( SET_LABEL(     ssCandsLabel, ps ) != "" );
  useBu   = ( SET_LABEL(     buCandsLabel, ps ) != "" );
  useBd   = ( SET_LABEL(     bdCandsLabel, ps ) != "" );
  useBs   = ( SET_LABEL(     bsCandsLabel, ps ) != "" );
  useK0   = ( SET_LABEL(     k0CandsLabel, ps ) != "" );
  useL0   = ( SET_LABEL(     l0CandsLabel, ps ) != "" );
  if ( useTrig ) consume< edm::TriggerResults           >( trigResultsToken,
                                                           trigResultsLabel );
  if ( useOnia ) consume< vector<pat::CompositeCandidate> >( oniaCandsToken,
                                                             oniaCandsLabel );
  if ( useSd   ) consume< vector<pat::CompositeCandidate> >(   sdCandsToken,
                                                               sdCandsLabel );
  if ( useSs   ) consume< vector<pat::CompositeCandidate> >(   ssCandsToken,
                                                               ssCandsLabel );
  if ( useBu   ) consume< vector<pat::CompositeCandidate> >(   buCandsToken,
                                                               buCandsLabel );
  if ( useBd   ) consume< vector<pat::CompositeCandidate> >(   bdCandsToken,
                                                               bdCandsLabel );
  if ( useBs   ) consume< vector<pat::CompositeCandidate> >(   bsCandsToken,
                                                               bsCandsLabel );
  if ( useK0   ) consume< vector<pat::CompositeCandidate> >(   k0CandsToken,
                                                               k0CandsLabel );
  if ( useL0   ) consume< vector<pat::CompositeCandidate> >(   l0CandsToken,
                                                               l0CandsLabel );

  static const BPHSoftMuonSelect sms;

  double  phiIMassMin =  0.85;
  double  phiIMassMax =  1.20;
  double  phiIPtMin   = 18.0;
  double  phiIEtaMax  = -1.0;
  double  phiIYMax    = -1.0;
  double  phiBMassMin =  0.85;
  double  phiBMassMax =  3.30;
  double  phiBPtMin   = 18.0;
  double  phiBEtaMax  = -1.0;
  double  phiBYMax    =  1.25;
  double jPsiIMassMin =  2.95;
  double jPsiIMassMax =  3.30;
  double jPsiIPtMin   = 25.0;
  double jPsiIEtaMax  = -1.0;
  double jPsiIYMax    = -1.0;
  double jPsiBMassMin =  2.95;
  double jPsiBMassMax =  3.30;
  double jPsiBPtMin   = 20.0;
  double jPsiBEtaMax  = -1.0;
  double jPsiBYMax    =  1.25;
  double psi2IMassMin =  3.40;
  double psi2IMassMax =  4.00;
  double psi2IPtMin   = 18.0;
  double psi2IEtaMax  = -1.0;
  double psi2IYMax    = -1.0;
  double psi2BMassMin =  3.40;
  double psi2BMassMax =  4.00;
  double psi2BPtMin   = 10.0;
  double psi2BEtaMax  = -1.0;
  double psi2BYMax    =  1.25;
  double  upsIMassMin =  8.50;
  double  upsIMassMax = 11.0;
  double  upsIPtMin   = 15.0;
  double  upsIEtaMax  = -1.0;
  double  upsIYMax    = -1.0;
  double  upsBMassMin =  8.50;
  double  upsBMassMax = 11.0;
  double  upsBPtMin   = 12.0;
  double  upsBEtaMax  =  1.5;
  double  upsBYMax    = -1.0;

  double oniaProbMin =  -0.005;
  double oniaCosMin  = -1.0;
  double oniaSigMin  = -1.0;

  double oniaMuPtMinLoose  =  2.0;
  double oniaMuPtMinTight  = -1.0;
  double oniaMuEtaMaxLoose = -1.0;
  double oniaMuEtaMaxTight = -1.0;

   phiIBasicSelect   = new BPHCompositeBasicSelect(
                           phiIMassMin,  phiIMassMax,
                           phiIPtMin  ,  phiIEtaMax ,  phiIYMax );
   phiBBasicSelect   = new BPHCompositeBasicSelect(
                           phiBMassMin,  phiBMassMax,
                           phiBPtMin  ,  phiBEtaMax ,  phiBYMax );
  jPsiIBasicSelect   = new BPHCompositeBasicSelect(
                           jPsiIMassMin, jPsiIMassMax,
                           jPsiIPtMin  , jPsiIEtaMax , jPsiIYMax );
  jPsiBBasicSelect   = new BPHCompositeBasicSelect(
                           jPsiBMassMin, jPsiBMassMax,
                           jPsiBPtMin  , jPsiBEtaMax , jPsiBYMax );
  psi2IBasicSelect   = new BPHCompositeBasicSelect(
                           psi2IMassMin, psi2IMassMax,
                           psi2IPtMin  , psi2IEtaMax , psi2IYMax );
  psi2BBasicSelect   = new BPHCompositeBasicSelect(
                           psi2BMassMin, psi2BMassMax,
                           psi2BPtMin  , psi2BEtaMax , psi2BYMax );
   upsIBasicSelect   = new BPHCompositeBasicSelect(
                           upsIMassMin,  upsIMassMax,
                           upsIPtMin  ,  upsIEtaMax ,  upsIYMax );
   upsBBasicSelect   = new BPHCompositeBasicSelect(
                           upsBMassMin,  upsBMassMax,
                           upsBPtMin  ,  upsBEtaMax ,  upsBYMax );
  oniaVertexSelect   = new BPHCompositeVertexSelect(
                           oniaProbMin, oniaCosMin, oniaSigMin );
  oniaDaughterSelect = new BPHDaughterSelect(
                           oniaMuPtMinLoose , oniaMuPtMinTight ,
                           oniaMuEtaMaxLoose, oniaMuEtaMaxTight, &sms );

  double npJPsiMassMin   = BPHParticleMasses::jPsiMass - 0.150;
  double npJPsiMassMax   = BPHParticleMasses::jPsiMass + 0.150;
  double npJPsiPtMin     =  8.0;
  double npJPsiEtaMax    = -1.0;
  double npJPsiYMax      = -1.0;
  double npMuPtMinLoose  =  4.0;
  double npMuPtMinTight  = -1.0;
  double npMuEtaMaxLoose =  2.2;
  double npMuEtaMaxTight = -1.0;

  npJPsiBasicSelect    = new BPHCompositeBasicSelect(
                             npJPsiMassMin, npJPsiMassMax,
                             npJPsiPtMin  , npJPsiEtaMax , npJPsiYMax );
  npJPsiDaughterSelect = new BPHDaughterSelect(
                             npMuPtMinLoose , npMuPtMinTight ,
                             npMuEtaMaxLoose, npMuEtaMaxTight, &sms );

  double buIMassMin       =      0.0;
  double buIMassMax       = 999999.0;
  double buIPtMin         = 27.0;
  double buIEtaMax        = -1.0;
  double buIYMax          = -1.0;
  double buIJPsiMassMin   = -1.0;//BPHParticleMasses::jPsiMass - 0.150;
  double buIJPsiMassMax   = -1.0;//BPHParticleMasses::jPsiMass + 0.150;
  double buIJPsiPtMin     = 25.0;
  double buIJPsiEtaMax    = -1.0;
  double buIJPsiYMax      = -1.0;
  double buIProbMin       =  0.15;
  double buICosMin        = -2.0;
  double buISigMin        = -1.0;
//  double buIMuPtMinLoose  = -1.0;
//  double buIMuPtMinTight  = -1.0;
//  double buIMuEtaMaxLoose = -1.0;
//  double buIMuEtaMaxTight = -1.0;

  buIKPtMin = 2.0;

  buIBasicSelect        = new BPHFittedBasicSelect(
                              buIMassMin, buIMassMax,
                              buIPtMin  , buIEtaMax , buIYMax );
  buIJPsiBasicSelect    = new BPHCompositeBasicSelect(
                              buIJPsiMassMin, buIJPsiMassMax,
                              buIJPsiPtMin  , buIJPsiEtaMax , buIJPsiYMax );
  buIVertexSelect       = new BPHVertexSelect( 'f',
                              buIProbMin, buICosMin, buISigMin );
  buIJPsiDaughterSelect = 0;
//  buIJPsiDaughterSelect = new BPHDaughterSelect(
//                              buIMuPtMinLoose , buIMuPtMinTight ,
//                              buIMuEtaMaxLoose, buMuEtaMaxTight, &sms );

  double buDMassMin       =      0.0;
  double buDMassMax       = 999999.0;
  double buDPtMin         = 10.0;
  double buDEtaMax        = -1.0;
  double buDYMax          = -1.0;
  double buDJPsiMassMin   = -1.0;//BPHParticleMasses::jPsiMass - 0.150;
  double buDJPsiMassMax   = -1.0;//BPHParticleMasses::jPsiMass + 0.150;
  double buDJPsiPtMin     =  8.0;
  double buDJPsiEtaMax    = -1.0;
  double buDJPsiYMax      = -1.0;
  double buDProbMin       =  0.10;
  double buDCosMin        =  0.99;
  double buDSigMin        =  3.0;
//  double buDMuPtMinLoose  = -1.0;
//  double buDMuPtMinTight  = -1.0;
//  double buDMuEtaMaxLoose = -1.0;
//  double buDMuEtaMaxTight = -1.0;

  buDKPtMin = 1.6;

  buDBasicSelect        = new BPHFittedBasicSelect(
                              buDMassMin, buDMassMax,
                              buDPtMin  , buDEtaMax , buDYMax );
  buDJPsiBasicSelect    = new BPHCompositeBasicSelect(
                              buDJPsiMassMin, buDJPsiMassMax,
                              buDJPsiPtMin  , buDJPsiEtaMax , buDJPsiYMax );
  buDVertexSelect       = new BPHVertexSelect( 'f',
                              buDProbMin, buDCosMin, buDSigMin );
  buDJPsiDaughterSelect = 0;
//  buDJPsiDaughterSelect = new BPHDaughterSelect(
//                              buDMuPtMinLoose , buDMuPtMinTight ,
//                              buDMuEtaMaxLoose, buDMuEtaMaxTight, &sms );

  double bdIMassMin       =      0.0;
  double bdIMassMax       = 999999.0;
  double bdIPtMin         = 27.0;
  double bdIEtaMax        = -1.0;
  double bdIYMax          = -1.0;
  double bdIJPsiMassMin   = -1.0;//BPHParticleMasses::jPsiMass - 0.150;
  double bdIJPsiMassMax   = -1.0;//BPHParticleMasses::jPsiMass + 0.150;
  double bdIJPsiPtMin     = 25.0;
  double bdIJPsiEtaMax    = -1.0;
  double bdIJPsiYMax      = -1.0;
  double bdIKx0MassMin    = BPHParticleMasses::kx0Mass - 0.100;
  double bdIKx0MassMax    = BPHParticleMasses::kx0Mass + 0.100;
  double bdIKx0PtMin      = -1.0;
  double bdIKx0EtaMax     = -1.0;
  double bdIKx0YMax       = -1.0;
  double bdIProbMin       =  0.10;
  double bdICosMin        =  0.99;
  double bdISigMin        =  3.0;
//  double bdIMuPtMinLoose  =  -1.0;
//  double bdIMuPtMinTight  =  -1.0;
//  double bdIMuEtaMaxLoose =  -1.0;
//  double bdIMuEtaMaxTight =  -1.0;

  bdIBasicSelect        = new BPHFittedBasicSelect(
                              bdIMassMin, bdIMassMax,
                              bdIPtMin  , bdIEtaMax , bdIYMax );
  bdIJPsiBasicSelect    = new BPHCompositeBasicSelect(
                              bdIJPsiMassMin, bdIJPsiMassMax,
                              bdIJPsiPtMin  , bdIJPsiEtaMax , bdIJPsiYMax );
  bdIKx0BasicSelect     = new BPHCompositeBasicSelect(
                              bdIKx0MassMin, bdIKx0MassMax,
                              bdIKx0PtMin  , bdIKx0EtaMax , bdIKx0YMax );
  bdIVertexSelect       = new BPHVertexSelect( 'f',
                              bdIProbMin, bdICosMin, bdISigMin );
  bdIJPsiDaughterSelect = 0;
//  bdIJPsiDaughterSelect = new BPHDaughterSelect(
//                              bdIMuPtMinLoose , bdIMuPtMinTight ,
//                              bdIMuEtaMaxLoose, bdIMuEtaMaxTight, &sms );

  double bdDMassMin       =      0.0;
  double bdDMassMax       = 999999.0;
  double bdDPtMin         = 10.0;
  double bdDEtaMax        = -1.0;
  double bdDYMax          = -1.0;
  double bdDJPsiMassMin   = -1.0;//BPHParticleMasses::jPsiMass - 0.150;
  double bdDJPsiMassMax   = -1.0;//BPHParticleMasses::jPsiMass + 0.150;
  double bdDJPsiPtMin     =  8.0;
  double bdDJPsiEtaMax    = -1.0;
  double bdDJPsiYMax      = -1.0;
  double bdDKx0MassMin    = BPHParticleMasses::kx0Mass - 0.100;
  double bdDKx0MassMax    = BPHParticleMasses::kx0Mass + 0.100;
  double bdDKx0PtMin      = -1.0;
  double bdDKx0EtaMax     = -1.0;
  double bdDKx0YMax       = -1.0;
  double bdDProbMin       =  0.10;
  double bdDCosMin        =  0.99;
  double bdDSigMin        =  3.0;
//  double bdDMuPtMinLoose  = -1.0;
//  double bdDMuPtMinTight  = -1.0;
//  double bdDMuEtaMaxLoose = -1.0;
//  double bdDMuEtaMaxTight = -1.0;

  bdDBasicSelect        = new BPHFittedBasicSelect(
                              bdDMassMin, bdDMassMax,
                              bdDPtMin  , bdDEtaMax , bdDYMax );
  bdDJPsiBasicSelect    = new BPHCompositeBasicSelect(
                              bdDJPsiMassMin, bdDJPsiMassMax,
                              bdDJPsiPtMin  , bdDJPsiEtaMax , bdDJPsiYMax );
  bdDKx0BasicSelect     = new BPHCompositeBasicSelect(
                              bdDKx0MassMin, bdDKx0MassMax,
                              bdDKx0PtMin  , bdDKx0EtaMax , bdDKx0YMax );
  bdDVertexSelect       = new BPHVertexSelect( 'f',
                              bdDProbMin, bdDCosMin, bdDSigMin );
  bdDJPsiDaughterSelect = 0;
//  bdDJPsiDaughterSelect = new BPHDaughterSelect(
//                              bdDMuPtMinLoose , bdDMuPtMinTight ,
//                              bdDMuEtaMaxLoose, bdDMuEtaMaxTight, &sms );

  double bsIMassMin       =      0.0;
  double bsIMassMax       = 999999.0;
  double bsIPtMin         = 27.0;
  double bsIEtaMax        = -1.0;
  double bsIYMax          = -1.0;
  double bsIJPsiMassMin   = BPHParticleMasses::jPsiMass - 0.150;
  double bsIJPsiMassMax   = BPHParticleMasses::jPsiMass + 0.150;
  double bsIJPsiPtMin     = 25.0;
  double bsIJPsiEtaMax    = -1.0;
  double bsIJPsiYMax      = -1.0;
  double bsIPhiMassMin    = BPHParticleMasses::phiMass - 0.010;
  double bsIPhiMassMax    = BPHParticleMasses::phiMass + 0.010;
  double bsIPhiPtMin      = -1.0;
  double bsIPhiEtaMax     = -1.0;
  double bsIPhiYMax       = -1.0;
  double bsIProbMin       =  0.10;
  double bsICosMin        =  0.99;
  double bsISigMin        =  3.0;
//  double bsIMuPtMinLoose  = -1.0;
//  double bsIMuPtMinTight  = -1.0;
//  double bsIMuEtaMaxLoose = -1.0;
//  double bsIMuEtaMaxTight = -1.0;

  bsIBasicSelect        = new BPHFittedBasicSelect(
                              bsIMassMin, bsIMassMax,
                              bsIPtMin  , bsIEtaMax , bsIYMax );
  bsIJPsiBasicSelect    = new BPHCompositeBasicSelect(
                              bsIJPsiMassMin, bsIJPsiMassMax,
                              bsIJPsiPtMin  , bsIJPsiEtaMax , bsIJPsiYMax );
  bsIPhiBasicSelect     = new BPHCompositeBasicSelect(
                              bsIPhiMassMin, bsIPhiMassMax,
                              bsIPhiPtMin  , bsIPhiEtaMax , bsIPhiYMax );
  bsIVertexSelect       = new BPHVertexSelect( 'f',
                              bsIProbMin, bsICosMin, bsISigMin );
  bsIJPsiDaughterSelect = 0;
//  bsIJPsiDaughterSelect = new BPHDaughterSelect(
//                              bsIMuPtMinLoose , bsIMuPtMinTight ,
//                              bsIMuEtaMaxLoose, bsIMuEtaMaxTight, &sms );

  double bsDMassMin       =      0.0;
  double bsDMassMax       = 999999.0;
  double bsDPtMin         = 10.0;
  double bsDEtaMax        = -1.0;
  double bsDYMax          = -1.0;
  double bsDJPsiMassMin   = -1.0;//BPHParticleMasses::jPsiMass - 0.150;
  double bsDJPsiMassMax   = -1.0;//BPHParticleMasses::jPsiMass + 0.150;
  double bsDJPsiPtMin     =  8.0;
  double bsDJPsiEtaMax    = -1.0;
  double bsDJPsiYMax      = -1.0;
  double bsDPhiMassMin    = BPHParticleMasses::phiMass - 0.010;
  double bsDPhiMassMax    = BPHParticleMasses::phiMass + 0.010;
  double bsDPhiPtMin      = -1.0;
  double bsDPhiEtaMax     = -1.0;
  double bsDPhiYMax       = -1.0;
  double bsDProbMin       =  0.10;
  double bsDCosMin        =  0.99;
  double bsDSigMin        =  3.0;
//  double bsDMuPtMinLoose  = -1.0;
//  double bsDMuPtMinTight  = -1.0;
//  double bsDMuEtaMaxLoose = -1.0;
//  double bsDMuEtaMaxTight = -1.0;

  bsDBasicSelect        = new BPHFittedBasicSelect(
                              bsDMassMin, bsDMassMax,
                              bsDPtMin  , bsDEtaMax , bsDYMax );
  bsDJPsiBasicSelect    = new BPHCompositeBasicSelect(
                              bsDJPsiMassMin, bsDJPsiMassMax,
                              bsDJPsiPtMin  , bsDJPsiEtaMax , bsDJPsiYMax );
  bsDPhiBasicSelect     = new BPHCompositeBasicSelect(
                              bsDPhiMassMin, bsDPhiMassMax,
                              bsDPhiPtMin  , bsDPhiEtaMax , bsDPhiYMax );
  bsDVertexSelect       = new BPHVertexSelect( 'f',
                              bsDProbMin, bsDCosMin, bsDSigMin );
  bsDJPsiDaughterSelect = 0;
//  bsDJPsiDaughterSelect = new BPHDaughterSelect(
//                              bsDMuPtMinLoose , bsDMuPtMinTight ,
//                              bsDMuEtaMaxLoose, bsDMuEtaMaxTight, &sms );

}


BPHHistoSpecificDecay::~BPHHistoSpecificDecay() {

  delete  phiIBasicSelect;
  delete  phiBBasicSelect;
  delete jPsiIBasicSelect;
  delete jPsiBBasicSelect;
  delete psi2IBasicSelect;
  delete psi2BBasicSelect;
  delete  upsIBasicSelect;
  delete  upsBBasicSelect;
  delete oniaVertexSelect;
  delete oniaDaughterSelect;

  delete npJPsiBasicSelect;
  delete npJPsiDaughterSelect;

  delete buIBasicSelect;
  delete buIJPsiBasicSelect;
  delete buIVertexSelect;
  delete buIJPsiDaughterSelect;
  delete buDBasicSelect;
  delete buDJPsiBasicSelect;
  delete buDVertexSelect;
  delete buDJPsiDaughterSelect;

  delete bdIBasicSelect;
  delete bdIJPsiBasicSelect;
  delete bdIKx0BasicSelect;
  delete bdIVertexSelect;
  delete bdIJPsiDaughterSelect;
  delete bdDBasicSelect;
  delete bdDJPsiBasicSelect;
  delete bdDKx0BasicSelect;
  delete bdDVertexSelect;
  delete bdDJPsiDaughterSelect;

  delete bsIBasicSelect;
  delete bsIJPsiBasicSelect;
  delete bsIPhiBasicSelect;
  delete bsIVertexSelect;
  delete bsIJPsiDaughterSelect;
  delete bsDBasicSelect;
  delete bsDJPsiBasicSelect;
  delete bsDPhiBasicSelect;
  delete bsDVertexSelect;
  delete bsDJPsiDaughterSelect;

}


void BPHHistoSpecificDecay::fillDescriptions(
                            edm::ConfigurationDescriptions& descriptions ) {
   edm::ParameterSetDescription desc;
   desc.add<string>( "trigResultsLabel", "" );
   desc.add<string>(   "oniaCandsLabel", "" );
   desc.add<string>(     "sdCandsLabel", "" );
   desc.add<string>(     "ssCandsLabel", "" );
   desc.add<string>(     "buCandsLabel", "" );
   desc.add<string>(     "bdCandsLabel", "" );
   desc.add<string>(     "bsCandsLabel", "" );
   desc.add<string>(     "k0CandsLabel", "" );
   desc.add<string>(     "l0CandsLabel", "" );
   descriptions.add( "process.bphHistoSpecificDecay", desc );
   return;
}


void BPHHistoSpecificDecay::beginJob() {

  createHisto( "massDIPhi"   ,  35, 0.85, 1.20 ); // Phi  mass inclusive
  createHisto( "massTIPhi"   ,  35, 0.85, 1.20 ); // Phi  mass inclusive
  createHisto( "massDBPhi"   ,  35, 0.85, 1.20 ); // Phi  mass barrel
  createHisto( "massTBPhi"   ,  35, 0.85, 1.20 ); // Phi  mass barrel
  createHisto( "massDIJPsi"  ,  35, 2.95, 3.30 ); // JPsi mass inclusive
  createHisto( "massTIJPsi"  ,  35, 2.95, 3.30 ); // JPsi mass inclusive
  createHisto( "massDBJPsi"  ,  35, 2.95, 3.30 ); // JPsi mass barrel
  createHisto( "massTBJPsi"  ,  35, 2.95, 3.30 ); // JPsi mass barrel
  createHisto( "massDIPsi2"  ,  60, 3.40, 4.00 ); // Psi2 mass inclusive
  createHisto( "massTIPsi2"  ,  60, 3.40, 4.00 ); // Psi2 mass inclusive
  createHisto( "massDBPsi2"  ,  60, 3.40, 4.00 ); // Psi2 mass barrel
  createHisto( "massTBPsi2"  ,  60, 3.40, 4.00 ); // Psi2 mass barrel
  createHisto( "massDIUps123", 125, 8.50, 11.0 ); // Ups  mass inclusive
  createHisto( "massTIUps123", 125, 8.50, 11.0 ); // Ups  mass inclusive
  createHisto( "massDBUps123", 125, 8.50, 11.0 ); // Ups  mass barrel
  createHisto( "massTBUps123", 125, 8.50, 11.0 ); // Ups  mass barrel
  createHisto( "massDIBu"    ,  50, 5.00, 6.00 ); // Bu   mass inclusive
  createHisto( "massTIBu"    ,  50, 5.00, 6.00 ); // Bu   mass inclusive
  createHisto( "massDDBu"    ,  50, 5.00, 6.00 ); // Bu   mass displaced
  createHisto( "massTDBu"    ,  50, 5.00, 6.00 ); // Bu   mass displaced
  createHisto( "massDIBd"    ,  50, 5.00, 6.00 ); // Bd   mass inclusive
  createHisto( "massTIBd"    ,  50, 5.00, 6.00 ); // Bd   mass inclusive
  createHisto( "massDDBd"    ,  50, 5.00, 6.00 ); // Bd   mass displaced
  createHisto( "massTDBd"    ,  50, 5.00, 6.00 ); // Bd   mass displaced
  createHisto( "massDIBs"    ,  50, 5.00, 6.00 ); // Bs   mass inclusive
  createHisto( "massTIBs"    ,  50, 5.00, 6.00 ); // Bs   mass inclusive
  createHisto( "massDDBs"    ,  50, 5.00, 6.00 ); // Bs   mass displaced
  createHisto( "massTDBs"    ,  50, 5.00, 6.00 ); // Bs   mass displaced
  createHisto( "mfitDIBu"    ,  50, 5.00, 6.00 ); // Bu   mass, with constraint
  createHisto( "mfitTIBu"    ,  50, 5.00, 6.00 ); // Bu   mass, with constraint
  createHisto( "mfitDDBu"    ,  50, 5.00, 6.00 ); // Bu   mass, with constraint
  createHisto( "mfitTDBu"    ,  50, 5.00, 6.00 ); // Bu   mass, with constraint
  createHisto( "mfitDIBd"    ,  50, 5.00, 6.00 ); // Bd   mass, with constraint
  createHisto( "mfitTIBd"    ,  50, 5.00, 6.00 ); // Bd   mass, with constraint
  createHisto( "mfitDDBd"    ,  50, 5.00, 6.00 ); // Bd   mass, with constraint
  createHisto( "mfitTDBd"    ,  50, 5.00, 6.00 ); // Bd   mass, with constraint
  createHisto( "mfitDIBs"    ,  50, 5.00, 6.00 ); // Bs   mass, with constraint
  createHisto( "mfitTIBs"    ,  50, 5.00, 6.00 ); // Bs   mass, with constraint
  createHisto( "mfitDDBs"    ,  50, 5.00, 6.00 ); // Bs   mass, with constraint
  createHisto( "mfitTDBs"    ,  50, 5.00, 6.00 ); // Bs   mass, with constraint
  createHisto( "massDBuJPsi" ,  35, 2.95, 3.30 ); // JPsi mass in Bu decay
  createHisto( "massTBuJPsi" ,  35, 2.95, 3.30 ); // JPsi mass in Bu decay
  createHisto( "massDBdJPsi" ,  35, 2.95, 3.30 ); // JPsi mass in Bd decay
  createHisto( "massTBdJPsi" ,  35, 2.95, 3.30 ); // JPsi mass in Bd decay
  createHisto( "massDBsJPsi" ,  35, 2.95, 3.30 ); // JPsi mass in Bs decay
  createHisto( "massTBsJPsi" ,  35, 2.95, 3.30 ); // JPsi mass in Bs decay
  createHisto( "massDBsPhi"  ,  50, 1.01, 1.03 ); // Phi  mass in Bs decay
  createHisto( "massTBsPhi"  ,  50, 1.01, 1.03 ); // Phi  mass in Bs decay
  createHisto( "massDBdKx0"  ,  50, 0.80, 1.05 ); // Kx0  mass in Bd decay
  createHisto( "massTBdKx0"  ,  50, 0.80, 1.05 ); // Kx0  mass in Bd decay
  createHisto( "massDK0s"    ,  50, 0.40, 0.60 ); // K0s  mass
  createHisto( "mfitDK0s"    ,  50, 0.40, 0.60 ); // K0s  mass
  createHisto( "massDLambda0",  50, 1.00, 1.30 ); // Lambda0 mass
  createHisto( "mfitDLambda0",  50, 1.00, 1.30 ); // Lambda0 mass

  createHisto( "massFull"   , 200, 2.00, 12.0 ); // Full onia mass

  createHisto( "ctauDIJPsi"  ,  60, -0.05, 0.25 ); // JPsi ctau inclusive
  createHisto( "ctauTIJPsi"  ,  60, -0.05, 0.25 ); // JPsi ctau inclusive
  createHisto( "ctauDBJPsi"  ,  60, -0.05, 0.25 ); // JPsi ctau barrel
  createHisto( "ctauTBJPsi"  ,  60, -0.05, 0.25 ); // JPsi ctau barrel
  createHisto( "ctauDIBu"    ,  60, -0.05, 0.25 ); // Bu   ctau inclusive
  createHisto( "ctauTIBu"    ,  60, -0.05, 0.25 ); // Bu   ctau inclusive
  createHisto( "ctauDDBu"    ,  60, -0.05, 0.25 ); // Bu   ctau displaced
  createHisto( "ctauTDBu"    ,  60, -0.05, 0.25 ); // Bu   ctau displaced
  createHisto( "ctauDIBd"    ,  60, -0.05, 0.25 ); // Bd   ctau inclusive
  createHisto( "ctauTIBd"    ,  60, -0.05, 0.25 ); // Bd   ctau inclusive
  createHisto( "ctauDDBd"    ,  60, -0.05, 0.25 ); // Bd   ctau displaced
  createHisto( "ctauTDBd"    ,  60, -0.05, 0.25 ); // Bd   ctau displaced
  createHisto( "ctauDIBs"    ,  60, -0.05, 0.25 ); // Bs   ctau inclusive
  createHisto( "ctauTIBs"    ,  60, -0.05, 0.25 ); // Bs   ctau inclusive
  createHisto( "ctauDDBs"    ,  60, -0.05, 0.25 ); // Bs   ctau displaced
  createHisto( "ctauTDBs"    ,  60, -0.05, 0.25 ); // Bs   ctau displaced

  return;

}

void BPHHistoSpecificDecay::analyze( const edm::Event& ev,
                                     const edm::EventSetup& es ) {

  // get magnetic field
  edm::ESHandle<MagneticField> magneticField;
  es.get<IdealMagneticFieldRecord>().get( magneticField );

  // get object collections
  // collections are got through "BPHTokenWrapper" interface to allow
  // uniform access in different CMSSW versions

  //////////// trigger results ////////////

  edm::Handle< edm::TriggerResults > trigResults;
  const edm::TriggerNames* trigNames = 0;
  if ( useTrig ) {
    trigResultsToken.get( ev, trigResults );
    if ( trigResults.isValid() )
         trigNames = &ev.triggerNames( *trigResults );
  }

  bool flag_Dimuon25_Jpsi                     = false;
  bool flag_Dimuon20_Jpsi_Barrel_Seagulls     = false;
  bool flag_Dimuon14_Phi_Barrel_Seagulls      = false;
  bool flag_Dimuon18_PsiPrime                 = false;
  bool flag_Dimuon10_PsiPrime_Barrel_Seagulls = false;
  bool flag_Dimuon12_Upsilon_eta1p5           = false;
  bool flag_DoubleMu4_JpsiTrk_Displaced       = false;
  if ( trigNames != 0 ) {
    const edm::TriggerNames::Strings& names = trigNames->triggerNames();
    int iObj;
    int nObj = names.size();
    for ( iObj = 0; iObj < nObj; ++iObj ) {
      CHK_TRIG( trigResults, names, iObj, Dimuon25_Jpsi                     )
      CHK_TRIG( trigResults, names, iObj, Dimuon20_Jpsi_Barrel_Seagulls     )
      CHK_TRIG( trigResults, names, iObj, Dimuon14_Phi_Barrel_Seagulls      )
      CHK_TRIG( trigResults, names, iObj, Dimuon18_PsiPrime                 )
      CHK_TRIG( trigResults, names, iObj, Dimuon10_PsiPrime_Barrel_Seagulls )
      CHK_TRIG( trigResults, names, iObj, Dimuon12_Upsilon_eta1p5           )
      CHK_TRIG( trigResults, names, iObj, DoubleMu4_JpsiTrk_Displaced       )
    }
  }

  //////////// quarkonia ////////////

  edm::Handle< vector<pat::CompositeCandidate> > oniaCands;
  int iqo;
  int nqo = 0;
  if ( useOnia ) {
    oniaCandsToken.get( ev, oniaCands );
    nqo = oniaCands->size();
  }

  for ( iqo = 0; iqo < nqo; ++ iqo ) {
    LogTrace( "DataDump" )
           << "*********** quarkonium " << iqo << "/" << nqo << " ***********";
    const pat::CompositeCandidate& cand = oniaCands->at( iqo );
//    if ( !oniaVertexSelect->accept( cand,
//                                    BPHUserData::getByRef<reco::Vertex>( cand,
//                                    "primaryVertex" ) ) ) continue;
    if ( !oniaDaughterSelect->accept( cand ) ) continue;
    fillHisto( "Full", cand, 'c' );
    if (  phiBBasicSelect->accept( cand ) ) {
      fillHisto( "DBPhi"   , cand, 'c' );
      if ( flag_Dimuon14_Phi_Barrel_Seagulls )
      fillHisto( "TBPhi"   , cand, 'c' );
    }
    if ( jPsiIBasicSelect->accept( cand ) ) {
      fillHisto( "DIJPsi"  , cand, 'c' );
      if ( flag_Dimuon25_Jpsi )
      fillHisto( "TIJPsi"  , cand, 'c' );
    }
    if ( jPsiBBasicSelect->accept( cand ) ) {
      fillHisto( "DBJPsi"  , cand, 'c' );
      if ( flag_Dimuon20_Jpsi_Barrel_Seagulls )
      fillHisto( "TBJPsi"  , cand, 'c' );
    }
    if ( psi2IBasicSelect->accept( cand ) ) {
      fillHisto( "DIPsi2"  , cand, 'c' );
      if ( flag_Dimuon18_PsiPrime )
      fillHisto( "TIPsi2"  , cand, 'c' );
    }
    if ( psi2BBasicSelect->accept( cand ) ) {
      fillHisto( "DBPsi2"  , cand, 'c' );
      if ( flag_Dimuon10_PsiPrime_Barrel_Seagulls )
      fillHisto( "TBPsi2"  , cand, 'c' );
    }
    if (  upsBBasicSelect->accept( cand ) ) {
      fillHisto( "DBUps123", cand, 'c' );
      if ( flag_Dimuon12_Upsilon_eta1p5 )
      fillHisto( "TBUps123", cand, 'c' );
    }
  }

  //////////// Bu ////////////

  edm::Handle< vector<pat::CompositeCandidate> > buCands;
  int ibu;
  int nbu = 0;
  if ( useBu ) {
    buCandsToken.get( ev, buCands );
    nbu = buCands->size();
  }

  for ( ibu = 0; ibu < nbu; ++ ibu ) {
    LogTrace( "DataDump" )
           << "*********** Bu " << ibu << "/" << nbu << " ***********";
    const pat::CompositeCandidate& cand = buCands->at( ibu );
    const pat::CompositeCandidate* jPsi = BPHUserData::getByRef
         <pat::CompositeCandidate>( cand, "refToJPsi" );
    LogTrace( "DataDump" )
           << "JPsi: " << jPsi;
    if ( jPsi == 0 ) continue;
    if ( !npJPsiBasicSelect   ->accept( *jPsi ) ) continue;
    if ( !npJPsiDaughterSelect->accept( *jPsi ) ) continue;
    const reco::Candidate* kptr = BPHDaughters::get( cand, 0.49, 0.50 ).front();
    if ( kptr == 0 ) continue;
    if ( buIBasicSelect       ->accept( *jPsi ) &&
         buIJPsiBasicSelect   ->accept( *jPsi ) &&
//         buIJPsiDaughterSelect->accept( *jPsi ) &&
         buIVertexSelect      ->accept( cand,
                                BPHUserData::getByRef<reco::Vertex>( *jPsi,
                                "primaryVertex" ) ) &&
         ( kptr->pt() > buIKPtMin ) ) {
      fillHisto( "DIBu"    ,  cand, 'f' );
      fillHisto( "DIBuJPsi", *jPsi, 'c' );
      if ( flag_Dimuon25_Jpsi ) {
	fillHisto( "TIBu"    ,  cand, 'f' );
	fillHisto( "TIBuJPsi", *jPsi, 'c' );
      }
    }
    if ( buDBasicSelect       ->accept( *jPsi ) &&
         buDJPsiBasicSelect   ->accept( *jPsi ) &&
//         buDJPsiDaughterSelect->accept( *jPsi ) &&
         buDVertexSelect      ->accept( cand,
                                BPHUserData::getByRef<reco::Vertex>( *jPsi,
                                "primaryVertex" ) ) &&
         ( kptr->pt() > buDKPtMin ) ) {
      fillHisto( "DDBu"    ,  cand, 'f' );
      fillHisto( "DDBuJPsi", *jPsi, 'c' );
      if ( flag_DoubleMu4_JpsiTrk_Displaced ) {
	fillHisto( "TDBu"    ,  cand, 'f' );
	fillHisto( "TDBuJPsi", *jPsi, 'c' );
      }
    }
  }

  //////////// Bd ////////////

  edm::Handle< vector<pat::CompositeCandidate> > bdCands;
  int ibd;
  int nbd = 0;
  if ( useBd ) {
    bdCandsToken.get( ev, bdCands );
    nbd = bdCands->size();
  }

  for ( ibd = 0; ibd < nbd; ++ ibd ) {
    LogTrace( "DataDump" )
           << "*********** Bd " << ibd << "/" << nbd << " ***********";
    const pat::CompositeCandidate& cand = bdCands->at( ibd );
    const pat::CompositeCandidate* jPsi = BPHUserData::getByRef
         <pat::CompositeCandidate>( cand, "refToJPsi" );
    LogTrace( "DataDump" )
           << "JPsi: " << jPsi;
    if ( jPsi == 0 ) continue;
    if ( !npJPsiBasicSelect   ->accept( *jPsi ) ) continue;
    if ( !npJPsiDaughterSelect->accept( *jPsi ) ) continue;
    const pat::CompositeCandidate* kx0 = BPHUserData::getByRef
         <pat::CompositeCandidate>( cand, "refToKx0" );
    LogTrace( "DataDump" )
           << "Kx0: " << kx0;
    if ( kx0 == 0 ) continue;
    if ( bdIBasicSelect       ->accept( *jPsi ) &&
         bdIJPsiBasicSelect   ->accept( *jPsi ) &&
         bdIKx0BasicSelect    ->accept( * kx0 ) &&
//         bdIJPsiDaughterSelect->accept( *jPsi ) &&
         bdIVertexSelect      ->accept( cand,
                                BPHUserData::getByRef<reco::Vertex>( *jPsi,
                                "primaryVertex" ) ) ) {
      fillHisto( "DIBd"    ,  cand, 'f' );
      fillHisto( "DIBdJPsi", *jPsi, 'c' );
      if ( flag_Dimuon25_Jpsi ) {
	fillHisto( "TIBd"    ,  cand, 'f' );
	fillHisto( "TIBdJPsi", *jPsi, 'c' );
      }
    }
    if ( bdDBasicSelect       ->accept( *jPsi ) &&
         bdDJPsiBasicSelect   ->accept( *jPsi ) &&
         bdDKx0BasicSelect    ->accept( * kx0 ) &&
//         bdDJPsiDaughterSelect->accept( *jPsi ) &&
         bdDVertexSelect      ->accept( cand,
                                BPHUserData::getByRef<reco::Vertex>( *jPsi,
                                "primaryVertex" ) ) ) {
      fillHisto( "DDBd"    ,  cand, 'f' );
      fillHisto( "DDBdJPsi", *jPsi, 'c' );
      if ( flag_DoubleMu4_JpsiTrk_Displaced ) {
	fillHisto( "TDBd"    ,  cand, 'f' );
	fillHisto( "TDBdJPsi", *jPsi, 'c' );
      }
    }
  }

  //////////// Bs ////////////

  edm::Handle< vector<pat::CompositeCandidate> > bsCands;
  int ibs;
  int nbs = 0;
  if ( useBs ) {
    bsCandsToken.get( ev, bsCands );
    nbs = bsCands->size();
  }

  for ( ibs = 0; ibs < nbs; ++ ibs ) {
    LogTrace( "DataDump" )
           << "*********** Bs " << ibs << "/" << nbs << " ***********";
    const pat::CompositeCandidate& cand = bsCands->at( ibs );
    const pat::CompositeCandidate* jPsi = BPHUserData::getByRef
         <pat::CompositeCandidate>( cand, "refToJPsi" );
    LogTrace( "DataDump" )
           << "JPsi: " << jPsi;
    if ( jPsi == 0 ) continue;
    const pat::CompositeCandidate* phi = BPHUserData::getByRef
         <pat::CompositeCandidate>( cand, "refToPhi" );
    LogTrace( "DataDump" )
           << "Phi: " << phi;
    if ( phi == 0 ) continue;
    if ( bsIBasicSelect       ->accept( *jPsi ) &&
         bsIJPsiBasicSelect   ->accept( *jPsi ) &&
         bsIPhiBasicSelect    ->accept( * phi ) &&
//         bsIJPsiDaughterSelect->accept( *jPsi ) &&
         bsIVertexSelect      ->accept( cand,
                                BPHUserData::getByRef<reco::Vertex>( *jPsi,
                                "primaryVertex" ) ) ) {
      fillHisto( "DIBs"    ,  cand, 'f' );
      fillHisto( "DIBsJPsi", *jPsi, 'c' );
      if ( flag_Dimuon25_Jpsi ) {
	fillHisto( "TIBs"    ,  cand, 'f' );
	fillHisto( "TIBsJPsi", *jPsi, 'c' );
      }
    }
    if ( bsDBasicSelect       ->accept( *jPsi ) &&
         bsDJPsiBasicSelect   ->accept( *jPsi ) &&
         bsDPhiBasicSelect    ->accept( * phi ) &&
//         bsDJPsiDaughterSelect->accept( *jPsi ) &&
         bsDVertexSelect      ->accept( cand,
                                BPHUserData::getByRef<reco::Vertex>( *jPsi,
                                "primaryVertex" ) ) ) {
      fillHisto( "DDBs"    ,  cand, 'f' );
      fillHisto( "DDBsJPsi", *jPsi, 'c' );
      if ( flag_DoubleMu4_JpsiTrk_Displaced ) {
	fillHisto( "TDBs"    ,  cand, 'f' );
	fillHisto( "TDBsJPsi", *jPsi, 'c' );
      }
    }
  }

  //////////// K0s ////////////

  edm::Handle< vector<pat::CompositeCandidate> > k0Cands;
  int ik0;
  int nk0 = 0;
  if ( useK0 ) {
    k0CandsToken.get( ev, k0Cands );
    nk0 = k0Cands->size();
  }

  for ( ik0 = 0; ik0 < nk0; ++ ik0 ) {
    LogTrace( "DataDump" )
           << "*********** K0 " << ik0 << "/" << nk0 << " ***********";
    const pat::CompositeCandidate& cand = k0Cands->at( ik0 );
    fillHisto( "DK0s", cand, 'f' );
  }

  //////////// Lambda0 ////////////

  edm::Handle< vector<pat::CompositeCandidate> > l0Cands;
  int il0;
  int nl0 = 0;
  if ( useL0 ) {
    l0CandsToken.get( ev, l0Cands );
    nl0 = l0Cands->size();
  }

  for ( il0 = 0; il0 < nl0; ++ il0 ) {
    LogTrace( "DataDump" )
           << "*********** Lambda0 " << il0 << "/" << nl0 << " ***********";
    const pat::CompositeCandidate& cand = l0Cands->at( il0 );
    fillHisto( "DLambda0", cand, 'f' );
  }

  return;

}


void BPHHistoSpecificDecay::endJob() {
  return;
}


void BPHHistoSpecificDecay::fillHisto( const string& name,
                                       const pat::CompositeCandidate& cand,
                                       char svType ) {
  float mass = ( cand.hasUserFloat( "fitMass" ) ?
                 cand.   userFloat( "fitMass" ) : -1 );
  fillHisto( "mass" + name, cand.mass() );
  fillHisto( "mfit" + name, mass );

  const reco::Vertex*
  pvtx = BPHUserData::getByRef<reco::Vertex>( cand, "primaryVertex" );
  if ( pvtx == 0 ) {
    const pat::CompositeCandidate* jPsi = BPHUserData::getByRef
         <pat::CompositeCandidate>( cand, "refToJPsi" );
    if ( jPsi == 0 ) return;
    pvtx = BPHUserData::getByRef<reco::Vertex>( *jPsi, "primaryVertex" );
  }
  if ( pvtx == 0 ) return;

  const reco::Vertex* svtx = 0;
  if ( svType == 'f' )
  svtx = BPHUserData::get<reco::Vertex>( cand, "fitVertex" );
  if ( svtx == 0 )
  svtx = BPHUserData::get<reco::Vertex>( cand, "vertex" );
  if ( svtx == 0 ) return;

  TVector3 disp2( svtx->x() - pvtx->x(),
                  svtx->y() - pvtx->y(),
                  0 );
//  TVector3 disp3( svtx->x() - pvtx->x(),
//                  svtx->y() - pvtx->y(),
//                  svtx->z() - pvtx->z() );
  TVector3 cmom2( cand.px(), cand.py(), 0 );
//  TVector3 cmom3( cand.px(), cand.py(), cand.pz() );
  float cosAlpha2 = disp2.Dot( cmom2 ) / ( disp2.Perp() * cmom2.Perp() );
//  float cosAlpha3 = disp3.Dot( cmom3 ) / ( disp3.Perp() * cmom3.Perp() );
  AlgebraicVector3 vmom2( cand.px(), cand.py(), 0 );
//  AlgebraicVector3 vmom3( cand.px(), cand.py(), cand.pz() );

  VertexDistanceXY vdistXY;
//  VertexDistance3D vdist3D;
  Measurement1D distXY = vdistXY.distance( *svtx, *pvtx );
//  Measurement1D dist3D = vdist3D.distance( *svtx, *pvtx );
//  float mass = cand.mass();
  double ctauPV2 = distXY.value() * cosAlpha2 * mass / cmom2.Perp();
/*
//  double ctauPV3 = dist3D.value() * cosAlpha3 * mass / cmom3.Perp();
  GlobalError sve = svtx->error();
  GlobalError pve = pvtx->error();
  AlgebraicSymMatrix33 dve = sve.matrix() + pve.matrix();
  double ctauErrPV2 = sqrt( ROOT::Math::Similarity( vmom2, dve ) ) * mass /
                            cmom2.Perp2();
//  double ctauErrPV3 = sqrt( ROOT::Math::Similarity( vmom3, dve ) ) * mass /
//                            cmom3.Perp2();
//  d2d = Measurement1D( ctauPV2, ctauPV2 / ctauErrPV2 );
//  d3d = Measurement1D( ctauPV3, ctauPV3 / ctauErrPV3 );
*/

  fillHisto( "ctau" + name, ctauPV2 );

  return;
}


void BPHHistoSpecificDecay::fillHisto( const string& name, float x ) {
  map<string,TH1F*>::iterator iter = histoMap.find( name );
  map<string,TH1F*>::iterator iend = histoMap.end();
  if ( iter == iend ) return;
  iter->second->Fill( x );
  return;
}


void BPHHistoSpecificDecay::createHisto( const string& name,
                                         int nbin, float hmin, float hmax ) {
  histoMap[name] = fs->make<TH1F>( name.c_str(), name.c_str(),
                                   nbin, hmin, hmax );
  return;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( BPHHistoSpecificDecay );
