#include "tb/data.hpp"

#include <algorithm>
#include <random>

namespace tb {

namespace {

constexpr std::size_t kBasesPerBlock = 32uz;

/* clang-format off */
constexpr std::uint8_t kNucleotideCoder[] = {
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255,   0, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255,   0,   1 ,  1,   0, 255, 255,   2,
      3, 255, 255,   2, 255,   1,   0, 255,
    255, 255,   0,   1,   3,   3,   2,   0,
    255,   3, 255, 255, 255, 255, 255, 255,
    255,   0,   1,   1,   0, 255, 255,   2,
      3, 255, 255,   2, 255,   1,   0, 255,
    255, 255,   0,   1,   3,   3,   2,   0,
    255,   3, 255, 255, 255, 255, 255, 255
};
/* clang-format on */

constexpr char kNucleotideDecoder[] = {'A', 'C', 'G', 'T'};

} // namespace

MockSequence::MockSequence(std::size_t n_bases, int seed)
    : n_bases_(n_bases),
      data_(n_bases / kBasesPerBlock + (n_bases < kBasesPerBlock)) {
  std::mt19937_64 rng_engine(seed);
  auto n_blocks = n_bases / kBasesPerBlock;
  std::uniform_int_distribution<std::uint64_t> distr;
  std::ranges::generate(data_,
                        [&] -> std::uint64_t { return distr(rng_engine); });
};

} // namespace tb
