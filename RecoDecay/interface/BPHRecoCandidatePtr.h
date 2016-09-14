#ifndef BPHAnalysis_RecoDecay_BPHRecoCandidatePtr_h
#define BPHAnalysis_RecoDecay_BPHRecoCandidatePtr_h

#include "BPHAnalysis/RecoDecay/interface/BPHGenericPtr.h"
class BPHRecoCandidate;
typedef BPHGenericPtr<      BPHRecoCandidate>::type BPHRecoCandidatePtr;
typedef BPHGenericPtr<const BPHRecoCandidate>::type BPHRecoConstCandPtr;

#endif
