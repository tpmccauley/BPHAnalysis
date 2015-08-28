#ifndef BPHDecayVertex_H
#define BPHDecayVertex_H
/** \class BPHDecayVertex
 *
 *  Description: 
 *     mid-level base class to reconstruct decay vertex
 *
 *  $Date: 2015-07-03 13:49:53 $
 *  $Revision: 1.1 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//----------------------
// Base Class Headers --
//----------------------
#include "BPHAnalysis/RecoDecay/interface/BPHDecayMomentum.h"

namespace edm {
  class EventSetup;
}

namespace reco {
  class TransientTrack;
  class Vertex;
}

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
#include "DataFormats/VertexReco/interface/Vertex.h"

//---------------
// C++ Headers --
//---------------
#include <vector>

//              ---------------------
//              -- Class Interface --
//              ---------------------

class BPHDecayVertex: public virtual BPHDecayMomentum {

 public:

  /** Constructor is protected
   *  this object can exist only as part of a derived class
   */

  /** Destructor
   */
  virtual ~BPHDecayVertex();

  /** Operations
   */

  /// check for valid reconstructed vertex
  virtual bool isValid() const;

  /// get reconstructed vertex
  virtual const reco::Vertex& vertex() const;

  /// get list of TransientTrack
  const std::vector<reco::TransientTrack>& transientTracks() const;

  /// get TransientTrack for a daughter
  reco::TransientTrack* getTransientTrack( const reco::Candidate* cand ) const;

 protected:

  // constructor
  BPHDecayVertex( const edm::EventSetup* es );

  // utility function used to cash reconstruction results
  virtual void setNotUpdated() const;

 private:

  // EventSetup needed to build TransientTrack
  const edm::EventSetup* evSetup;

  // reconstruction results cache
  mutable bool oldTracks;
  mutable bool oldVertex;
  mutable bool validVertex;
  mutable std::vector<reco::TransientTrack> trTracks;
  mutable std::map<const reco::Candidate*,reco::TransientTrack*> ttMap;
  mutable reco::Vertex fittedVertex;

  // create TransientTrack and fit vertex
  virtual void tTracks() const;
  virtual void fitVertex() const;

};


#endif // BPHDecayVertex_H

