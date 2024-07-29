#pragma once

#include "tb/data.hpp"

namespace tb {

struct MinimizeArgs {
  MockSequence const& seq;
  std::int32_t window_length;
  std::int32_t kmer_length;
};

// Reference naive implementations
std::vector<KMer> NaiveMinimize(MinimizeArgs);

// Deque based implementations
std::vector<KMer> DequeMinimize(MinimizeArgs);
std::vector<KMer> InplaceMinimize(MinimizeArgs);
std::vector<KMer> RingMinimize(MinimizeArgs);

// Arg min based implementations
std::vector<KMer> ArgMinMinimize(MinimizeArgs);
std::vector<KMer> ArgMinEveMinimize(MinimizeArgs);
std::vector<KMer> ArgMinDuffMinimize(MinimizeArgs);
std::vector<KMer> ArgMinUnrolledMinimize(MinimizeArgs);
std::vector<KMer> NtHashArgMinMinimize(MinimizeArgs);
std::vector<KMer> NtHashPrecomputedArgMinMinimize(MinimizeArgs);
std::vector<KMer> NtHashPrecomputedArgMinUnrolledMinimize(MinimizeArgs);

// Arg min recovery based implementations
std::vector<KMer> ArgMinRecoveryMinimize(MinimizeArgs);
std::vector<KMer> ArgMinRecoveryEveMinimize(MinimizeArgs);
std::vector<KMer> ArgMinRecoveryDuffMinimize(MinimizeArgs);
std::vector<KMer> ArgMinRecoveryUnrolledMinimize(MinimizeArgs);
std::vector<KMer> NtHashArgMinRecoveryMinimize(MinimizeArgs);
std::vector<KMer> NtHashPrecomputedArgMinRecoveryMinimize(MinimizeArgs);
std::vector<KMer> NtHashPrecomputedArgMinUnrolledRecoveryMinimize(MinimizeArgs);

// Arg min rolling
std::vector<KMer> ArgMinRollingMinimize(MinimizeArgs);

}  // namespace tb
