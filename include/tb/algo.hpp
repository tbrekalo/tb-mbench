#pragma once

#include "tb/data.hpp"

namespace tb {

struct MinimizeArgs {
  MockSequence const& seq;
  std::int32_t window_length;
  std::int32_t kmer_length;
};

std::vector<KMer> NaiveMinimize(MinimizeArgs);
std::vector<KMer> DequeMinimize(MinimizeArgs);
std::vector<KMer> InplaceMinimize(MinimizeArgs);
std::vector<KMer> RingMinimize(MinimizeArgs);
std::vector<KMer> ArgMinMinimize(MinimizeArgs);
std::vector<KMer> NtHashArgMinMinimize(MinimizeArgs);
std::vector<KMer> ArgMinRecoveryMinimize(MinimizeArgs);
std::vector<KMer> ArgMinRecoveryEveMinimize(MinimizeArgs);
std::vector<KMer> NtHashArgMinRecoveryMinimize(MinimizeArgs);

}  // namespace tb
