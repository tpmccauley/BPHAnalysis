#ifndef BPHAnalysis_SpecificDecay_BPHBcToJPsiPiBuilder_h
#define BPHAnalysis_SpecificDecay_BPHBcToJPsiPiBuilder_h
/** \class BPHBcToJPsiPiBuilder
 *
 *  Description: 
 *     Class to build Bc to JPsi pi+- candidates
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

class BPHParticleNeutralVeto;
class BPHParticlePtSelect;
class BPHParticleEtaSelect;
class BPHMassSelect;
class BPHChi2Select;
class BPHMassFitSelect;

//---------------
// C++ Headers --
//---------------
#include <string>
#include <vector>

//              ---------------------
//              -- Class Interface --
//              ---------------------

class BPHBcToJPsiPiBuilder {

 public:

  /** Constructor
   */
  BPHBcToJPsiPiBuilder( const edm::EventSetup& es,
      const std::vector<BPHPlusMinusConstCandPtr>& jpsiCollection,
      const BPHRecoBuilder::BPHGenericCollection*  pionCollection );

  /** Destructor
   */
  virtual ~BPHBcToJPsiPiBuilder();

  /** Operations
   */
  /// build Bu candidates
  std::vector<BPHRecoConstCandPtr> build();

  /// set cuts
  void setPiPtMin    ( double pt  );
  void setPiEtaMax   ( double eta );
  void setJPsiMassMin( double m   );
  void setJPsiMassMax( double m   );
  void setMassMin    ( double m   );
  void setMassMax    ( double m   );
  void setProbMin    ( double p   );
  void setMassFitMin ( double m   );
  void setMassFitMax ( double m   );
  void setConstr     ( bool flag  );

  /// get current cuts
  double getPiPtMin    () const;
  double getPiEtaMax   () const;
  double getJPsiMassMin() const;
  double getJPsiMassMax() const;
  double getMassMin    () const;
  double getMassMax    () const;
  double getProbMin    () const;
  double getMassFitMin () const;
  double getMassFitMax () const;
  bool   getConstr     () const;

 private:

  // private copy and assigment constructors
  BPHBcToJPsiPiBuilder           ( const BPHBcToJPsiPiBuilder& x );
  BPHBcToJPsiPiBuilder& operator=( const BPHBcToJPsiPiBuilder& x );

  std::string jPsiName;
  std::string pionName;

  const edm::EventSetup* evSetup;
  const std::vector<BPHPlusMinusConstCandPtr>* jCollection;
  const BPHRecoBuilder::BPHGenericCollection*  pCollection;

  BPHMassSelect         * jpsiSel;
  BPHParticleNeutralVeto* pnVeto;
  BPHParticlePtSelect   *   ptSel;
  BPHParticleEtaSelect  *  etaSel;

  BPHMassSelect         * massSel;
  BPHChi2Select         * chi2Sel;
  BPHMassFitSelect      * mFitSel;

  bool massConstr;
  float minPDiff;
  bool updated;

  std::vector<BPHRecoConstCandPtr> bcList;

};


#endif

