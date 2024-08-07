#pragma once

#include "tb/data.hpp"

namespace tb {

struct MinimizeArgs {
  MockSequence const& seq;
  std::int32_t window_length;
  std::int32_t kmer_length;
};

std::vector<KMer::value_type> NtHash(MinimizeArgs);
std::vector<KMer::value_type> NtHashOpt(MinimizeArgs);

// Reference naive implementations
std::vector<KMer> NaiveMinimize(MinimizeArgs);

// Deque based implementations
std::vector<KMer> DequeMinimize(MinimizeArgs);

// Arg min based implementations
std::vector<KMer> ArgMinMinimize(MinimizeArgs);
std::vector<KMer> ArgMinUnrolledMinimize(MinimizeArgs);
std::vector<KMer> NtHashArgMinUnrolledMinimize(MinimizeArgs);

// Arg min recovery based implementations
std::vector<KMer> ArgMinRecoveryMinimize(MinimizeArgs);
std::vector<KMer> ArgMinRecoveryUnrolledMinimize(MinimizeArgs);
std::vector<KMer> NtHashRecoveryUnrolledMinimize(MinimizeArgs);

// Split window
std::vector<KMer> SplitWindowMinimize(MinimizeArgs);

}  // namespace tb
