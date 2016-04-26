#ifndef BPHPlusMinusConstCandPtr_H
#define BPHPlusMinusConstCandPtr_H

#include <memory>
class BPHPlusMinusCandidate;
typedef std::shared_ptr<      BPHPlusMinusCandidate> BPHPlusMinusCandidatePtr;
typedef std::shared_ptr<const BPHPlusMinusCandidate> BPHPlusMinusConstCandPtr;

#endif
