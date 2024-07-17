#pragma once

#include "tb/data.hpp"

namespace tb {

struct MinimizeArgs {
  MockSequence const &seq;
  std::int32_t window_length;
  std::int32_t kmer_length;
};

std::vector<KMer> NaiveMinimize(MinimizeArgs);
std::vector<KMer> DequeMinimize(MinimizeArgs);
std::vector<KMer> InplaceMinimize(MinimizeArgs);
std::vector<KMer> RingMinimize(MinimizeArgs);
std::vector<KMer> ArgminMinimize(MinimizeArgs);
std::vector<KMer> ArgminEveMinimize(MinimizeArgs);

} // namespace tb
