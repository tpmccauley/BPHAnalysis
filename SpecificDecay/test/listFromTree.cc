#include <string>
#include <iostream>
#include <fstream>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

using namespace std;

class BPHReco {

 public :

  TTree* currentTree;

  UInt_t   runNumber;
  UInt_t   lumiSection;
  UInt_t   eventNumber;
  string*  recoName;
  Float_t  recoMass;
  Float_t  recoTime;
  Float_t  recoErrT;

  TBranch* b_runNumber;
  TBranch* b_lumiSection;
  TBranch* b_eventNumber;
  TBranch* b_recoName;
  TBranch* b_recoMass;
  TBranch* b_recoTime;
  TBranch* b_recoErrT;

  class OutList {
   public:
    ofstream* osPtr;
    bool printMass;
    bool printTime;
    bool printErrT;
    OutList( ofstream* op, bool mass = true,
                           bool time = false,
                           bool errt = false ):
      osPtr( op ), printMass( mass ), printTime( time ), printErrT( errt ) {
    }
  };
  map<string,OutList*> ofMap;

  BPHReco( const string& head, const string& tail ) {
    recoName = new string( "" );
    ofMap["TBPhi"    ] = new OutList( new ofstream( ( head + "BarPhi"
                                                           + tail ).c_str() ) );
    ofMap["TIJPsi"   ] = new OutList( new ofstream( ( head + "IncJPsi"
                                                           + tail ).c_str() ) );
    ofMap["TBJPsi"   ] = new OutList( new ofstream( ( head + "BarJPsi"
                                                           + tail ).c_str() ) );
    ofMap["TIPsi2"   ] = new OutList( new ofstream( ( head + "IncPsi2"
                                                           + tail ).c_str() ) );
    ofMap["TBPsi2"   ] = new OutList( new ofstream( ( head + "BarPsi2"
                                                           + tail ).c_str() ) );
    ofMap["TBUps123" ] = new OutList( new ofstream( ( head + "BarUpsilon123"
                                                           + tail ).c_str() ) );
    ofMap["TIBu"     ] = new OutList( new ofstream( ( head + "InclusiveBu"
                                                           + tail ).c_str() ),
                                      true, true, true );
    ofMap["TDBu"     ] = new OutList( new ofstream( ( head + "DisplacedBu"
                                                           + tail ).c_str() ) );
//    ofMap["TIBd"     ] = new OutList( new ofstream( ( head + "InclusiveBd"
//                                                           + tail ).c_str() ) );
    ofMap["TDBd"     ] = new OutList( new ofstream( ( head + "DisplacedBd"
                                                           + tail ).c_str() ) );
//    ofMap["TIBs"     ] = new OutList( new ofstream( ( head + "InclusiveBs"
//                                                           + tail ).c_str() ) );
    ofMap["TDBs"     ] = new OutList( new ofstream( ( head + "DisplacedBs"
                                                           + tail ).c_str() ) );
//    ofMap["DK0s"     ] = new OutList( new ofstream( ( head + "K0s"
//                                                           + tail ).c_str() ) );
//    ofMap["DLambda0" ] = new OutList( new ofstream( ( head + "Lambda0"
//                                                           + tail ).c_str() ) );
    ofMap["TIB0"     ] = new OutList( new ofstream( ( head + "InclusiveB0"
                                                           + tail ).c_str() ) );
//    ofMap["TDB0"     ] = new OutList( new ofstream( ( head + "DisplacedB0"
//                                                           + tail ).c_str() ) );
    ofMap["TILambdab"] = new OutList( new ofstream( ( head + "InclusiveLambdab"
                                                           + tail ).c_str() ) );
    ofMap["TDLambdab"] = new OutList( new ofstream( ( head + "DisplacedLambdab"
                                                           + tail ).c_str() ) );
    ofMap["TIBc"     ] = new OutList( new ofstream( ( head + "InclusiveBc"
                                                           + tail ).c_str() ) );
    ofMap["TDBc"     ] = new OutList( new ofstream( ( head + "DisplacedBc"
                                                           + tail ).c_str() ) );
    ofMap["TIX3872"  ] = new OutList( new ofstream( ( head + "InclusiveX3872"
                                                           + tail ).c_str() ) );
    ofMap["TDX3872"  ] = new OutList( new ofstream( ( head + "DisplacedX3872"
                                                           + tail ).c_str() ) );
  }
  ~BPHReco() {}

  void init( TTree* tree ) {
    currentTree = tree;
    tree->SetMakeClass( 1 );
    tree->SetBranchAddress( "runNumber"  , &runNumber  , &b_runNumber   );
    tree->SetBranchAddress( "lumiSection", &lumiSection, &b_lumiSection );
    tree->SetBranchAddress( "eventNumber", &eventNumber, &b_eventNumber );
    tree->SetBranchAddress( "recoName"   , &recoName   , &b_recoName );
    tree->SetBranchAddress( "recoMass"   , &recoMass   , &b_recoMass );
    tree->SetBranchAddress( "recoTime"   , &recoTime   , &b_recoTime );
    tree->SetBranchAddress( "recoErrT"   , &recoErrT   , &b_recoErrT );
    return;
  }

  void process() {
    map<string,OutList*>::const_iterator iter = ofMap.find( *recoName );
    map<string,OutList*>::const_iterator iend = ofMap.end();
    ofstream* osPtr = 0;
    if ( iter != iend ) {
      const OutList* ol = iter->second;
      if ( ( recoErrT > 0.0 ) || !ol->printErrT ) {
        osPtr = ol->osPtr;
        if ( ol->printMass ) *osPtr << ' ' << recoMass;
        if ( ol->printTime ) *osPtr << ' ' << recoTime;
        if ( ol->printErrT ) *osPtr << ' ' << recoErrT;
        *osPtr << endl;
      }
    }
    return;
  }

  void loop() {
    // get number of events
    int evtnum = static_cast<int>( currentTree->GetEntries() );
    int ientry;
    for ( ientry = 0; ientry < evtnum; ientry++ ) {
      currentTree->GetEntry( ientry );
//      if ( runNumber < 315975 ) continue;
      process();
    }
    return;
  }

  void loop( const string& treeListName ) {
    ifstream treeListFile( treeListName.c_str() );
    char* treeListLine = new char[1000];
    char* treeListLptr;
    char* treeFileName;
    while ( treeListFile.getline( treeListLine, 1000 ) ) {
      treeListLptr = treeListLine;
      while ( *treeListLptr == ' ' ) treeListLptr++;
      treeFileName = treeListLptr;
      while ( ( *treeListLptr != ' ' ) &&
              ( *treeListLptr != '\0' ) ) treeListLptr++;
      *treeListLptr = '\0';
      cout << "open file " << treeFileName << endl;
      TChain* c = new TChain( "bphHistoSpecificDecay/BPHReco" );
      c->Add( treeFileName );
      init( c );
      loop();
      delete c;
    }
  }

};

int main( int argc, char* argv[] ) {
  string list( argv[1] );
  string head( argv[2] );
  string tail( argv[3] );
  BPHReco r( head, tail );
  r.loop( list );
  return 0;
}

