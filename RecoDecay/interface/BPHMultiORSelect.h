#ifndef BPHMultiORSelect_H
#define BPHMultiORSelect_H
/** \class BPHMultiORSelect
 *
 *  Description: 
 *     Class to combine multiple selection (OR mode)
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
class BPHMultiORSelect: public T {

 public:

  /** Constructor
   */
  BPHMultiORSelect();

  /** Destructor
   */
  virtual ~BPHMultiORSelect();

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
  BPHMultiORSelect           ( const BPHMultiORSelect<T>& x );
  BPHMultiORSelect& operator=( const BPHMultiORSelect<T>& x );

  struct SelectElement {
    T* selector;
    bool mode;
  };

  std::vector<SelectElement> selectList;

};

#include "BPHAnalysis/RecoDecay/interface/BPHMultiORSelect.hpp"

#endif // BPHMultiORSelect_H

