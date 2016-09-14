#ifndef BPHAnalysis_RecoDecay_BPHPlusMinusCandidatePtr_h
#define BPHAnalysis_RecoDecay_BPHPlusMinusCandidatePtr_h

#include "BPHAnalysis/RecoDecay/interface/BPHGenericPtr.h"
class BPHPlusMinusCandidate;
typedef BPHGenericPtr<      BPHPlusMinusCandidate>::type
                            BPHPlusMinusCandidatePtr;
typedef BPHGenericPtr<const BPHPlusMinusCandidate>::type
                            BPHPlusMinusConstCandPtr;

#endif
