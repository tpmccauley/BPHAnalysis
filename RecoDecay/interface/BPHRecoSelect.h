#ifndef BPHRecoSelect_H
#define BPHRecoSelect_H
/** \class BPHRecoSelect
 *
 *  Description: 
 *     Base class for daughter particle selection
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
namespace reco {
  class Candidate;
}

//---------------
// C++ Headers --
//---------------
#include <string>

//              ---------------------
//              -- Class Interface --
//              ---------------------

class BPHRecoSelect {

 public:

  /** Constructor
   */
  BPHRecoSelect();

  /** Destructor
   */
  virtual ~BPHRecoSelect();

  /** Operations
   */
  /// accept function
  /// pointers to other particles in the decays can be obtained 
  /// by the function "get" giving the particle name
  virtual bool accept( const reco::Candidate& cand ) const = 0;

 protected:

  // function to get other particles pointers
  const reco::Candidate* get( const std::string& name ) const;

 private:

  // private copy and assigment constructors
  BPHRecoSelect           ( const BPHRecoSelect& x );
  BPHRecoSelect& operator=( const BPHRecoSelect& x );

};


#endif // BPHRecoSelect_H

