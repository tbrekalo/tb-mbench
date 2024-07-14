#pragma once

#include "tb/data.hpp"

namespace tb {

struct MinimizeArgs {
  MockSequence const &seq;
  std::uint32_t kmer_length;
  std::uint32_t window_length;
};

std::vector<KMer> NaiveMinimize(MinimizeArgs);

} // namespace tb
