#ifndef TestBaseNtuple_h
#define TestBaseNtuple_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include <string>

class BPHRecoCandidate;

class TestBPHRecoDecay: public edm::EDAnalyzer {

 public:

  explicit TestBPHRecoDecay( const edm::ParameterSet& ps );
  virtual ~TestBPHRecoDecay();

  virtual void beginJob();
  virtual void analyze( const edm::Event& ev, const edm::EventSetup& es );
  virtual void endJob();

 private:

  void dumpRecoCand( const std::string& name,
                     const BPHRecoCandidate* cand );

};

#endif
