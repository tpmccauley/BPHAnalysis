#ifndef BPHAnalysis_RecoDecay_BPHGenericPtr_h
#define BPHAnalysis_RecoDecay_BPHGenericPtr_h

#include <memory>
template <class T>
class BPHGenericPtr {
 public:
  typedef typename std::shared_ptr<T> type;
};

#endif
