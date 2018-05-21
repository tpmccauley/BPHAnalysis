#ifndef BPHAnalysis_SpecificDecay_BPHBdToJPsiKsBuilder_h
#define BPHAnalysis_SpecificDecay_BPHBdToJPsiKsBuilder_h
/** \class BPHBdToJPsiKsBuilder
 *
 *  Description: 
 *     Class to build B0 to JPsi K0s candidates
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

class BPHBdToJPsiKsBuilder {

 public:

  /** Constructor
   */
  BPHBdToJPsiKsBuilder( const edm::EventSetup& es,
      const std::vector<BPHPlusMinusConstCandPtr>& jpsiCollection,
      const std::vector<BPHPlusMinusConstCandPtr>&  k0sCollection );

  /** Destructor
   */
  virtual ~BPHBdToJPsiKsBuilder();

  /** Operations
   */
  /// build Bs candidates
  std::vector<BPHRecoConstCandPtr> build();

  /// get original daughters map
  const std::map<const BPHRecoCandidate*,
                 const BPHRecoCandidate*>& daughMap() const;

  /// set cuts
  void setJPsiMassMin( double m  );
  void setJPsiMassMax( double m  );
  void setK0MassMin  ( double m  );
  void setK0MassMax  ( double m  );
  void setMassMin    ( double m  );
  void setMassMax    ( double m  );
  void setProbMin    ( double p  );
  void setMassFitMin ( double m  );
  void setMassFitMax ( double m  );
  void setConstr     ( bool flag );

  /// get current cuts
  double getJPsiMassMin() const;
  double getJPsiMassMax() const;
  double getK0MassMin  () const;
  double getK0MassMax  () const;
  double getMassMin    () const;
  double getMassMax    () const;
  double getProbMin    () const;
  double getMassFitMin () const;
  double getMassFitMax () const;
  bool   getConstr     () const;

 private:

  // private copy and assigment constructors
  BPHBdToJPsiKsBuilder           ( const BPHBdToJPsiKsBuilder& x );
  BPHBdToJPsiKsBuilder& operator=( const BPHBdToJPsiKsBuilder& x );

  std::string jPsiName;
  std::string  k0sName;

  const edm::EventSetup* evSetup;
  const std::vector<BPHPlusMinusConstCandPtr>* jCollection;
  const std::vector<BPHPlusMinusConstCandPtr>* kCollection;

  BPHMassSelect   * jpsiSel;
  BPHMassFitSelect* mk0sSel;

  BPHMassSelect   * massSel;
  double probMin;
  BPHMassFitSelect* mFitSel;

  bool massConstr;
  float minPDiff;
  bool updated;

  std::map<const BPHRecoCandidate*,const BPHRecoCandidate*> dMap;
  std::vector<BPHRecoConstCandPtr> bdList;

};


#endif

