#ifndef BPHAnalysis_SpecificDecay_BPHPhiToKKBuilder_h
#define BPHAnalysis_SpecificDecay_BPHPhiToKKBuilder_h
/** \class BPHPhiToKKBuilder
 *
 *  Description: 
 *     Class to build Phi to K+ K- candidates
 *
 *  \author Paolo Ronchese INFN Padova
 *
 */

//----------------------
// Base Class Headers --
//----------------------


//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
#include "BPHAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "BPHAnalysis/RecoDecay/interface/BPHRecoCandidate.h"
#include "BPHAnalysis/RecoDecay/interface/BPHPlusMinusCandidate.h"

#include "FWCore/Framework/interface/Event.h"

class BPHParticlePtSelect;
class BPHParticleEtaSelect;
class BPHChi2Select;
class BPHMassSelect;

//---------------
// C++ Headers --
//---------------
#include <string>
#include <vector>

//              ---------------------
//              -- Class Interface --
//              ---------------------

class BPHPhiToKKBuilder {

 public:

  /** Constructor
   */
  BPHPhiToKKBuilder( const edm::EventSetup& es,
       const BPHRecoBuilder::BPHGenericCollection* pCollection,
       const BPHRecoBuilder::BPHGenericCollection* nCollection );

  /** Destructor
   */
  virtual ~BPHPhiToKKBuilder();

  /** Operations
   */
  /// build Phi candidates
  std::vector<BPHPlusMinusConstCandPtr> build();

  /// set cuts
  void setPtMin  ( double pt  );
  void setEtaMax ( double eta );
  void setMassMin( double m   );
  void setMassMax( double m   );
  void setProbMin( double p   );
// mass constraint: foreseen but actually not implemented
  void setConstr ( double mass, double sigma );

  /// get current cuts
  double getPtMin  () const;
  double getEtaMax () const;
  double getMassMin() const;
  double getMassMax() const;
  double getProbMin() const;
// mass constraint: foreseen but actually not implemented
  double getConstrMass () const;
  double getConstrSigma() const;

 private:

  // private copy and assigment constructors
  BPHPhiToKKBuilder           ( const BPHPhiToKKBuilder& x );
  BPHPhiToKKBuilder& operator=( const BPHPhiToKKBuilder& x );

  std::string kPosName;
  std::string kNegName;

  const edm::EventSetup* evSetup;
  const BPHRecoBuilder::BPHGenericCollection* posCollection;
  const BPHRecoBuilder::BPHGenericCollection* negCollection;

  double ptMin;
  double etaMax;
  double massMin;
  double massMax;
  double probMin;

// old code left for example/reference

//  BPHParticlePtSelect *  ptSel;
//  BPHParticleEtaSelect* etaSel;
//  BPHMassSelect* massSel;
//  BPHChi2Select* chi2Sel;

//  class DZSelect;
//  DZSelect* dzSel;

// mass constraint: foreseen but actually not implemented
  double cMass;
  double cSigma;

  bool updated;

  std::vector<BPHPlusMinusConstCandPtr> phiList;

};


#endif

