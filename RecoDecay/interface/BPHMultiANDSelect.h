#ifndef BPHMultiANDSelect_H
#define BPHMultiANDSelect_H
/** \class BPHMultiANDSelect
 *
 *  Description: 
 *     Class to combine multiple selection (AND mode)
 *
 *
 *  $Date: 2015-07-06 18:40:19 $
 *  $Revision: 1.1 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//----------------------
// Base Class Headers --
//----------------------


//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class BPHRecoBuilder;
class BPHDecayMomentum;
class BPHDecayVertex;

namespace reco {
  class Candidate;
}

//---------------
// C++ Headers --
//---------------
#include <vector>

//              ---------------------
//              -- Class Interface --
//              ---------------------

template<class T>
class BPHMultiANDSelect: public T {

 public:

  /** Constructor
   */
  BPHMultiANDSelect();

  /** Destructor
   */
  virtual ~BPHMultiANDSelect();

  /** Operations
   */
  /// include selection
  void include( T& s, bool m = true );

  /// accept function
  virtual bool accept( const reco::Candidate & cand,
                       const BPHRecoBuilder*  build ) const;
  virtual bool accept( const reco::Candidate & cand ) const;
  virtual bool accept( const BPHDecayMomentum& cand ) const;
  virtual bool accept( const BPHDecayVertex  & cand ) const;

 private:

  // private copy and assigment constructors
  BPHMultiANDSelect           ( const BPHMultiANDSelect<T>& x );
  BPHMultiANDSelect& operator=( const BPHMultiANDSelect<T>& x );

  struct SelectElement {
    T* selector;
    bool mode;
  };

  std::vector<SelectElement> selectList;

};

#include "BPHAnalysis/RecoDecay/interface/BPHMultiANDSelect.hpp"

#endif // BPHMultiANDSelect_H

