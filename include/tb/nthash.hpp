#pragma once
#include <array>
#include <cstdint>
#include <limits>

#include "tb/data.hpp"

namespace tb {

// 64-bit random seeds corresponding to bases and their complements
const std::uint64_t kSeedA = 0x3c8b'fbb3'95c6'0470;
const std::uint64_t kSeedC = 0x3193'c185'62a0'2b4c;
const std::uint64_t kSeedG = 0x2032'3ed0'8257'2324;
const std::uint64_t kSeedT = 0x2d2a'04e6'7531'0c18;

inline constexpr std::array<std::uint64_t, 4> kNtHashSeeds = {
    kSeedA,
    kSeedC,
    kSeedG,
    kSeedT,
};

constexpr auto srol = []<class... Args> [[using gnu: always_inline, const]] (
                          Args&&... args) -> KMer::value_type {
  auto single_impl = [](std::uint64_t x) noexcept {
    /* clang-format off */
    uint64_t m = ((x & 0x8000000000000000ULL) >> 30) |
                 ((x & 0x100000000ULL) >> 32);
    /* clang-format on */
    return ((x << 1) & 0xFFFFFFFDFFFFFFFFULL) | m;
  };

  auto multi_impl = [] [[using gnu: always_inline]] (
                        KMer::value_type x, KMer::value_type n_rotations) {
    if (n_rotations == 0) {
      return x;
    }
    uint64_t v = (x << n_rotations) | (x >> (64 - n_rotations));
    uint64_t y = (v ^ (v >> 33)) &
                 (std::numeric_limits<uint64_t>::max() >> (64 - n_rotations));
    return v ^ (y | (y << 33));
  };

  if constexpr (requires(Args... args) { single_impl(args...); }) {
    return single_impl(std::forward<Args>(args)...);
  } else if constexpr (requires(Args... args) { multi_impl(args...); }) {
    return multi_impl(std::forward<Args>(args)...);
  } else {
    static_assert(false, "Overload not supported");
  }
};

constexpr auto nthash = [] [[using gnu: const, hot]]
                        (KMer::value_type prev, std::uint8_t base_out,
                         std::uint8_t base_in, std::uint64_t k) {
                          return srol(prev) ^ srol(kNtHashSeeds[base_out], k) ^
                                 kNtHashSeeds[base_in];
                        };

}  // namespace tb
