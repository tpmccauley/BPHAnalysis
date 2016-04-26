#ifndef AnalyzerTokenWrapper_H
#define AnalyzerTokenWrapper_H
/** \class TokenWrapper
 *
 *  Description: 
 *    common interface to get objects from "old" and "new" CMSSW version
 *    in an uniform way
 *
 *  $Date: 2016-04-15 17:47:56 $
 *  $Revision: 1.1 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//----------------------
// Base Class Headers --
//----------------------
#include "FWCore/Framework/interface/EDAnalyzer.h"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
#include "FWCore/Framework/interface/Event.h"

//---------------
// C++ Headers --
//---------------
#include <string>

//              ---------------------
//              -- Class Interface --
//              ---------------------

template<class Obj>
class TokenWrapper {
 public:
  typedef typename edm::InputTag type;
  bool get( const edm::Event& ev,
            edm::Handle<Obj>& obj ) {
    return ev.getByLabel( token, obj );
  }
  type token;
};

class AnalyzerWrapper: public edm::EDAnalyzer {
 protected:
  template<class Obj>
  void consume( TokenWrapper<Obj>& token,
                const std::string& label ) {
    token.token = edm::InputTag( label );
    return;
  }
};

#endif

