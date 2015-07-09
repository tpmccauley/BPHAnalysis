/*
 *  See header file for a description of this class.
 *
 *  $Date: 2015-07-06 18:40:19 $
 *  $Revision: 1.1 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//-----------------------
// This Class' Header --
//-----------------------
#include "BPHAnalysis/RecoDecay/interface/BPHRecoSelect.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BPHAnalysis/RecoDecay/interface/BPHRecoBuilder.h"

//---------------
// C++ Headers --
//---------------
#include <map>

using namespace std;

//-------------------
// Initializations --
//-------------------


//----------------
// Constructors --
//----------------
BPHRecoSelect::BPHRecoSelect() {
}

//--------------
// Destructor --
//--------------
BPHRecoSelect::~BPHRecoSelect() {
}

//--------------
// Operations --
//--------------
const reco::Candidate* BPHRecoSelect::get( const std::string& name ) const {
  map<std::string,
      const reco::Candidate*>& cMap = BPHRecoBuilder::daugMap();
  map<std::string,
      const reco::Candidate*>::iterator iter = cMap.find( name );
  return ( iter != cMap.end() ? iter->second : 0 );
}

