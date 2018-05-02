#ifndef BPHAnalysis_SpecificDecay_BPHHistoSpecificDecay_h
#define BPHAnalysis_SpecificDecay_BPHHistoSpecificDecay_h

#include "BPHAnalysis/RecoDecay/interface/BPHAnalyzerTokenWrapper.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Ref.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <string>

class TH1F;
class TVector3;

namespace reco {
  class Candidate;
  class Vertex;
}

class BPHHistoSpecificDecay:
      public BPHAnalyzerWrapper<BPHModuleWrapper::one_analyzer> {

 public:

  explicit BPHHistoSpecificDecay( const edm::ParameterSet& ps );
  virtual ~BPHHistoSpecificDecay();

  static void fillDescriptions( edm::ConfigurationDescriptions& descriptions );

  virtual void beginJob();
  virtual void analyze( const edm::Event& ev, const edm::EventSetup& es );
  virtual void endJob();

  class CandidateSelect {
   public:
    virtual ~CandidateSelect() {}
    virtual bool accept( const pat::CompositeCandidate& cand,
                         const reco::Vertex* pv = 0 ) const = 0 ;
  };

 private:

  std::string trigResultsLabel;
  std::string   oniaCandsLabel;
  std::string     sdCandsLabel;
  std::string     ssCandsLabel;
  std::string     buCandsLabel;
  std::string     bdCandsLabel;
  std::string     bsCandsLabel;
  std::string     k0CandsLabel;
  std::string     l0CandsLabel;
  std::string     b0CandsLabel;
  std::string     lbCandsLabel;
  BPHTokenWrapper< edm::TriggerResults > trigResultsToken;
  BPHTokenWrapper< std::vector<pat::CompositeCandidate> > oniaCandsToken;
  BPHTokenWrapper< std::vector<pat::CompositeCandidate> >   sdCandsToken;
  BPHTokenWrapper< std::vector<pat::CompositeCandidate> >   ssCandsToken;
  BPHTokenWrapper< std::vector<pat::CompositeCandidate> >   buCandsToken;
  BPHTokenWrapper< std::vector<pat::CompositeCandidate> >   bdCandsToken;
  BPHTokenWrapper< std::vector<pat::CompositeCandidate> >   bsCandsToken;
  BPHTokenWrapper< std::vector<pat::CompositeCandidate> >   k0CandsToken;
  BPHTokenWrapper< std::vector<pat::CompositeCandidate> >   l0CandsToken;
  BPHTokenWrapper< std::vector<pat::CompositeCandidate> >   b0CandsToken;
  BPHTokenWrapper< std::vector<pat::CompositeCandidate> >   lbCandsToken;
  bool useTrig;
  bool useOnia;
  bool useSd;
  bool useSs;
  bool useBu;
  bool useBd;
  bool useBs;
  bool useK0;
  bool useL0;
  bool useB0;
  bool useLb;

  edm::Service<TFileService> fs;
  std::map<std::string,TH1F*> histoMap;

  CandidateSelect*  phiIBasicSelect;
  CandidateSelect* jPsiIBasicSelect;
  CandidateSelect* psi2IBasicSelect;
  CandidateSelect*  upsIBasicSelect;
  CandidateSelect*  phiBBasicSelect;
  CandidateSelect* jPsiBBasicSelect;
  CandidateSelect* psi2BBasicSelect;
  CandidateSelect*  upsBBasicSelect;
  CandidateSelect* oniaVertexSelect;
  CandidateSelect* oniaDaughterSelect;

  CandidateSelect* npJPsiBasicSelect;
  CandidateSelect* npJPsiDaughterSelect;

  CandidateSelect* buIBasicSelect;
  CandidateSelect* buIJPsiBasicSelect;
  CandidateSelect* buIVertexSelect;
  CandidateSelect* buIJPsiDaughterSelect;
  CandidateSelect* buDBasicSelect;
  CandidateSelect* buDJPsiBasicSelect;
  CandidateSelect* buDVertexSelect;
  CandidateSelect* buDJPsiDaughterSelect;

  CandidateSelect* bdIBasicSelect;
  CandidateSelect* bdIJPsiBasicSelect;
  CandidateSelect* bdIKx0BasicSelect;
  CandidateSelect* bdIVertexSelect;
  CandidateSelect* bdIJPsiDaughterSelect;
  CandidateSelect* bdDBasicSelect;
  CandidateSelect* bdDJPsiBasicSelect;
  CandidateSelect* bdDKx0BasicSelect;
  CandidateSelect* bdDVertexSelect;
  CandidateSelect* bdDJPsiDaughterSelect;

  CandidateSelect* bsIBasicSelect;
  CandidateSelect* bsIJPsiBasicSelect;
  CandidateSelect* bsIPhiBasicSelect;
  CandidateSelect* bsIVertexSelect;
  CandidateSelect* bsIJPsiDaughterSelect;
  CandidateSelect* bsDBasicSelect;
  CandidateSelect* bsDJPsiBasicSelect;
  CandidateSelect* bsDPhiBasicSelect;
  CandidateSelect* bsDVertexSelect;
  CandidateSelect* bsDJPsiDaughterSelect;

  CandidateSelect* b0IBasicSelect;
  CandidateSelect* b0IJPsiBasicSelect;
  CandidateSelect* b0IK0sBasicSelect;
  CandidateSelect* b0IVertexSelect;
  CandidateSelect* b0IJPsiDaughterSelect;
  CandidateSelect* b0DBasicSelect;
  CandidateSelect* b0DJPsiBasicSelect;
  CandidateSelect* b0DK0sBasicSelect;
  CandidateSelect* b0DVertexSelect;
  CandidateSelect* b0DJPsiDaughterSelect;

  CandidateSelect* lbIBasicSelect;
  CandidateSelect* lbIJPsiBasicSelect;
  CandidateSelect* lbILambda0BasicSelect;
  CandidateSelect* lbIVertexSelect;
  CandidateSelect* lbIJPsiDaughterSelect;
  CandidateSelect* lbDBasicSelect;
  CandidateSelect* lbDJPsiBasicSelect;
  CandidateSelect* lbDLambda0BasicSelect;
  CandidateSelect* lbDVertexSelect;
  CandidateSelect* lbDJPsiDaughterSelect;

  double buIKPtMin;
  double buDKPtMin;

  void fillHisto   ( const std::string& name,
                     const pat::CompositeCandidate& cand, char svType );
  void fillHisto   ( const std::string& name, float x );
  void createHisto ( const std::string& name,
                     int nbin, float hmin, float hmax );

};

#endif
