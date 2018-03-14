#ifndef BPHAnalysis_SpecificDecay_BPHMassSelect_h
#define BPHAnalysis_SpecificDecay_BPHMassSelect_h
/** \class BPHMassSelect
 *
 *  Description: 
 *     Class for candidate selection by invariant mass (at momentum sum level)
 *
 *  \author Paolo Ronchese INFN Padova
 *
 */

//----------------------
// Base Class Headers --
//----------------------
#include "BPHAnalysis/RecoDecay/interface/BPHMomentumSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHMassCuts.h"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
#include "BPHAnalysis/RecoDecay/interface/BPHDecayMomentum.h"

//---------------
// C++ Headers --
//---------------


//              ---------------------
//              -- Class Interface --
//              ---------------------

class BPHMassSelect: public BPHMomentumSelect, public BPHMassCuts {

 public:

  /** Constructor
   */
  BPHMassSelect( double minMass, double maxMass ): BPHMassCuts( minMass,
                                                                maxMass ) {}

  /** Destructor
   */
  virtual ~BPHMassSelect() {}

  /** Operations
   */
  /// select particle
  virtual bool accept( const BPHDecayMomentum& cand ) const {
    double mass = cand.composite().mass();
    return ( ( mass >= mMin ) && ( mass <= mMax ) );
  }

 private:

  // private copy and assigment constructors
  BPHMassSelect           ( const BPHMassSelect& x );
  BPHMassSelect& operator=( const BPHMassSelect& x );

};


#endif

