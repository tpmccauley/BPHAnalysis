#ifndef BPHKinematicFit_H
#define BPHKinematicFit_H
/** \class BPHKinematicFit
 *
 *  Description: 
 *
 *
 *  $Date: 2015-07-03 17:02:36 $
 *  $Revision: 1.1 $
 *  \author Paolo Ronchese INFN Padova
 *     high-level base class to perform a kinematic fit
 */

//----------------------
// Base Class Headers --
//----------------------
#include "BPHAnalysis/RecoDecay/interface/BPHDecayVertex.h"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"

//---------------
// C++ Headers --
//---------------
#include <map>

//              ---------------------
//              -- Class Interface --
//              ---------------------

class BPHKinematicFit: public virtual BPHDecayVertex {

 public:

  /** Constructor is protected
   *  this object can exist only as part of a derived class
   */

  /** Destructor
   */
  virtual ~BPHKinematicFit();

  /** Operations
   */
  /// 
  /// add a simple particle giving it a name
  /// particles are cloned, eventually specifying a different mass
  /// and a sigma
  virtual void add( const std::string& name,
                    const reco::Candidate* daug,
                    double mass = -1.0, double sigma = -1.0 );

  /// apply a mass constraint
  void setConstraint( double mass, double sigma );
  /// retrieve the constraint
  double constrMass() const;
  double constrSigma() const;

  /// perform the kinematic fit and get the result
  virtual const RefCountedKinematicTree& kinematicTree() const;

  /// compute total momentum after the fit
  virtual const math::XYZTLorentzVector& p4() const;

 protected:

  // constructors
  BPHKinematicFit();
  // pointer used to retrieve informations from other bases
  BPHKinematicFit( const BPHKinematicFit* ptr );

  // utility function used to cash reconstruction results
  virtual void setNotUpdated() const;

 private:

  // mass constraint
  double massConst;
  double massSigma;

  // map linking daughters to mass sigma
  std::map<const reco::Candidate*,double> dMSig;

  // reconstruction results cache
  mutable bool updatedFit;
  mutable bool updatedMom;
  mutable RefCountedKinematicTree* kinTree;
  mutable math::XYZTLorentzVector totalMomentum;

  // per form the fit and compute the total momentum
  virtual void kinFit() const;
  virtual void fitMomentum() const;

};


#endif // BPHKinematicFit_H

