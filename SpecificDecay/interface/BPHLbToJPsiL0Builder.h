#ifndef BPHAnalysis_SpecificDecay_BPHLbToJPsiL0Builder_h
#define BPHAnalysis_SpecificDecay_BPHLbToJPsiL0Builder_h
/** \class BPHLbToJPsiL0Builder
 *
 *  Description: 
 *     Class to build Lambda_b to JPsi Lambda_0 candidates
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

class BPHLbToJPsiL0Builder {

 public:

  /** Constructor
   */
  BPHLbToJPsiL0Builder( const edm::EventSetup& es,
      const std::vector<BPHPlusMinusConstCandPtr>& jpsiCollection,
      const std::vector<BPHPlusMinusConstCandPtr>&   l0Collection );

  /** Destructor
   */
  virtual ~BPHLbToJPsiL0Builder();

  /** Operations
   */
  /// build Bs candidates
  std::vector<BPHRecoConstCandPtr> build();

  /// set cuts
  void setJPsiMassMin   ( double m  );
  void setJPsiMassMax   ( double m  );
  void setLambda0MassMin( double m  );
  void setLambda0MassMax( double m  );
  void setMassMin       ( double m  );
  void setMassMax       ( double m  );
  void setProbMin       ( double p  );
  void setMassFitMin    ( double m  );
  void setMassFitMax    ( double m  );
  void setConstr        ( bool flag );

  /// get current cuts
  double getJPsiMassMin   () const;
  double getJPsiMassMax   () const;
  double getLambda0MassMin() const;
  double getLambda0MassMax() const;
  double getMassMin       () const;
  double getMassMax       () const;
  double getProbMin       () const;
  double getMassFitMin    () const;
  double getMassFitMax    () const;
  bool   getConstr        () const;

 private:

  // private copy and assigment constructors
  BPHLbToJPsiL0Builder           ( const BPHLbToJPsiL0Builder& x );
  BPHLbToJPsiL0Builder& operator=( const BPHLbToJPsiL0Builder& x );

  std::string jPsiName;
  std::string   l0Name;

  const edm::EventSetup* evSetup;
  const std::vector<BPHPlusMinusConstCandPtr>* jCollection;
  const std::vector<BPHPlusMinusConstCandPtr>* lCollection;

  BPHMassSelect   * jpsiSel;
  BPHMassSelect   *  ml0Sel;

  BPHMassSelect   * massSel;
  double probMin;
  BPHMassFitSelect* mFitSel;

  bool massConstr;
  float minPDiff;
  bool updated;

  std::vector<BPHRecoConstCandPtr> lbList;

};


#endif

