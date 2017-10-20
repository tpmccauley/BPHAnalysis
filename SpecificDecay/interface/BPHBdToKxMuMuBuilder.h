#ifndef BPHAnalysis_SpecificDecay_BPHBdToKxMuMuBuilder_h
#define BPHAnalysis_SpecificDecay_BPHBdToKxMuMuBuilder_h
/** \class BPHBdToKxMuMuBuilder
 *
 *  Description: 
 *     Class to build B0 to K*0 mu+ mu- candidates
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

class BPHBdToKxMuMuBuilder {

 public:

  /** Constructor
   */
  BPHBdToKxMuMuBuilder( const edm::EventSetup& es,
      const std::vector<BPHPlusMinusConstCandPtr>& oniaCollection,
      const std::vector<BPHPlusMinusConstCandPtr>&  kx0Collection );

  /** Destructor
   */
  virtual ~BPHBdToKxMuMuBuilder();

  /** Operations
   */
  /// build Bs candidates
  std::vector<BPHRecoConstCandPtr> build();

  /// set cuts
  void setOniaMassMin( double m  );
  void setOniaMassMax( double m  );
  void setKxMassMin  ( double m  );
  void setKxMassMax  ( double m  );
  void setMassMin    ( double m  );
  void setMassMax    ( double m  );
  void setProbMin    ( double p  );
  void setMassFitMin ( double m  );
  void setMassFitMax ( double m  );
  void setConstr     ( bool flag );

  /// get current cuts
  double getOniaMassMin() const;
  double getOniaMassMax() const;
  double getKxMassMin  () const;
  double getKxMassMax  () const;
  double getMassMin    () const;
  double getMassMax    () const;
  double getProbMin    () const;
  double getMassFitMin () const;
  double getMassFitMax () const;
  bool   getConstr     () const;

 private:

  // private copy and assigment constructors
  BPHBdToKxMuMuBuilder           ( const BPHBdToKxMuMuBuilder& x );
  BPHBdToKxMuMuBuilder& operator=( const BPHBdToKxMuMuBuilder& x );

  std::string oniaName;
  std::string  kx0Name;

  const edm::EventSetup* evSetup;
  const std::vector<BPHPlusMinusConstCandPtr>* jCollection;
  const std::vector<BPHPlusMinusConstCandPtr>* kCollection;

  BPHMassSelect   * oniaSel;
  BPHMassSelect   * mkx0Sel;

  BPHMassSelect   * massSel;
  BPHChi2Select   * chi2Sel;
  BPHMassFitSelect* mFitSel;

  bool massConstr;
  float minPDiff;
  bool updated;

  std::vector<BPHRecoConstCandPtr> bdList;

};


#endif

