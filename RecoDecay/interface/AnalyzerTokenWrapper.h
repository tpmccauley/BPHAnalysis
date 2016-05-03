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
#include "FWCore/Framework/interface/EDConsumerBase.h"

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
  typedef typename edm::EDGetTokenT<Obj> type;
  bool get( const edm::Event& ev,
            edm::Handle<Obj>& obj ) {
    return ev.getByToken( token, obj );
  }
  type token;
};

template<class T>
class AnalyzerWrapper: public T {
 protected:
  template<class Obj>
  void consume( TokenWrapper<Obj>& token,
                const std::string& label ) {
    edm::InputTag tag( label );
    token.token = this->template consumes<Obj>( tag );
    return;
  }
};

#endif

