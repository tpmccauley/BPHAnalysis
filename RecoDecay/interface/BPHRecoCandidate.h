#ifndef BPHRecoCandidate_H
#define BPHRecoCandidate_H
/** \class BPHRecoCandidate
 *
 *  Description: 
 *     High level class for reconstructed decay candidates:
 *     - only constructor interfaces are defined in this class
 *     - functions returning results are defined in base classes:
 *       BPHDecayMomentum : decay products simple momentum sum
 *       BPHDecayVertex   : vertex reconstruction
 *       BPHKinematicFit  : kinematic fit and fitted momentum sum
 *
 *  $Date: 2015-07-03 10:08:22 $
 *  $Revision: 1.1 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//----------------------
// Base Class Headers --
//----------------------
#include "BPHAnalysis/RecoDecay/interface/BPHKinematicFit.h"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
#include "BPHAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

namespace edm {
  class EventSetup;
}

namespace reco {
  class Candidate;
}

//---------------
// C++ Headers --
//---------------
#include <vector>
#include <set>

//              ---------------------
//              -- Class Interface --
//              ---------------------

class BPHRecoCandidate: public virtual BPHKinematicFit {

 public:

  /** Constructors are private
   */
  /// create an "empty" object to add daughters later
  /// (see BPHDecayMomentum)
  BPHRecoCandidate( const edm::EventSetup* es );
  // create an object with daughters as specified in the ComponentSet
  BPHRecoCandidate( const edm::EventSetup* es,
                    const BPHRecoBuilder::ComponentSet& compSet );

  /** Destructor
   */
  virtual ~BPHRecoCandidate();

  /** Operations
   */
  /// look for candidates starting from particle collections as
  /// specified in the BPHRecoBuilder
  static std::vector<const BPHRecoCandidate*> build(
                                              const BPHRecoBuilder& builder,
                                              double mass = -1,
                                              double msig = -1 );

  /// delete all the candidates
  /// (all candidates are created in the heap)
  static void clear();

 protected:

  // template function called by "build" to allow
  // the creation of derived objects
  template <class T>
  static void fill( std::vector<const T*>& cList,
                    const BPHRecoBuilder& builder,
                    double mass = -1, double msig = -1 );

 private:

  // candidate bookkeeping to allow cleanup
  static std::set<const BPHRecoCandidate*> allCand;
  void addCand();

};


template <class T>
void BPHRecoCandidate::fill( std::vector<const T*>& cList,
                             const BPHRecoBuilder& builder,
                             double mass, double msig ) {
  // create paricle combinations
  const std::vector<BPHRecoBuilder::ComponentSet> dll = builder.build();
  // loop over combinations and create reconstructed particles
  int i;
  int n = dll.size();
  cList.reserve( n );
  T* rc = 0;
  for ( i = 0; i < n; ++i ) {
    // create reconstructed particle
    rc = new T( builder.eventSetup(), dll[i] );
    // apply mass constraint, if requested
    if ( mass > 0 ) rc->setConstraint( mass, msig );
    // apply post selection
    if ( builder.accept( *rc ) ) cList.push_back( rc );
    else                         delete rc;
  }
  return;
}

#endif // BPHRecoCandidate_H

