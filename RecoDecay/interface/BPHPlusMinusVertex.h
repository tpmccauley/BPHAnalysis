#ifndef BPHPlusMinusVertex_H
#define BPHPlusMinusVertex_H
/** \class BPHPlusMinusVertex
 *
 *  Description: 
 *     class for reconstructed decay vertices to opposite charge
 *     particle pairs
 *
 *  $Date: 2015-07-03 17:45:25 $
 *  $Revision: 1.1 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//----------------------
// Base Class Headers --
//----------------------
#include "BPHAnalysis/RecoDecay/interface/BPHDecayVertex.h"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

//---------------
// C++ Headers --
//---------------
#include <iostream>

//              ---------------------
//              -- Class Interface --
//              ---------------------

class BPHPlusMinusVertex: public virtual BPHDecayVertex {

 public:

  /** Constructor is protected
   *  this object can exist only as part of a derived class
   */

  /** Destructor
   */
  virtual ~BPHPlusMinusVertex();

  /** Operations
   */
  /// compute distance of closest approach
  virtual ClosestApproachInRPhi& cAppInRPhi() const;

 protected:

  BPHPlusMinusVertex( const edm::EventSetup* es );

  // utility functions to check/enforce the number of decay particles
  // at 2
  template<class T> static
  bool chkName( const T& cont,
                const std::string& name,
                const std::string& msg );
  template<class T> static
  bool chkSize( const T& cont,
                const std::string& msg );
  bool chkSize( const std::string& msg ) const;

  // utility function used to cash reconstruction results
  virtual void setNotUpdated() const;

 private:

  // reconstruction results cache
  mutable bool updated;
  mutable ClosestApproachInRPhi* inRPhi;

  // compute closest approach distance and cache it
  virtual void computeApp() const;

};


template<class T>
bool BPHPlusMinusVertex::chkName( const T& cont,
                                  const std::string& name,
                                  const std::string& msg ) {
  if ( cont.find( name ) != cont.end() ) return true;
  std::cout << msg << ", " << name << " not found" << std::endl;
  return false;
}


template<class T>
bool BPHPlusMinusVertex::chkSize( const T& cont,
                                  const std::string& msg ) {
  int n = cont.size();
  if ( n == 2 ) return true;
  std::cout << msg << ", size = " << n << std::endl;
  return false;
}

#endif // BPHPlusMinusVertex_H

