#pragma once

#include "tb/data.hpp"

namespace tb {

struct MinimizeArgs {
  MockSequence const &seq;
  std::uint32_t window_length;
  std::uint32_t kmer_length;
};

std::vector<KMer> NaiveMinimize(MinimizeArgs);
std::vector<KMer> DequeMinimize(MinimizeArgs);
std::vector<KMer> InplaceMinimize(MinimizeArgs);
std::vector<KMer> RingMinimize(MinimizeArgs);

} // namespace tb
