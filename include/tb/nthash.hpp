#pragma once
#include <immintrin.h>

#include <array>
#include <cstdint>
#include <limits>
#include <utility>

#include "immintrin.h"
#include "tb/data.hpp"

namespace tb {

// 64-bit random seeds corresponding to bases and their complements
inline constexpr std::uint64_t kSeedA = 0x3c8b'fbb3'95c6'0470;
inline constexpr std::uint64_t kSeedC = 0x3193'c185'62a0'2b4c;
inline constexpr std::uint64_t kSeedG = 0x2032'3ed0'8257'2324;
inline constexpr std::uint64_t kSeedT = 0x2d2a'04e6'7531'0c18;

inline constexpr std::array<std::uint64_t, 4> kNtHashSeeds = {
    kSeedA,
    kSeedC,
    kSeedG,
    kSeedT,
};

inline constexpr auto srol =
    []<class... Args> [[using gnu: always_inline, const]] (
        Args&&... args) -> KMer::value_type {
  auto single_impl =
      [] [[using gnu: always_inline, const]] (std::uint64_t x) noexcept {
        /* clang-format off */
    uint64_t m = ((x & 0x8000000000000000ULL) >> 30) |
                 ((x & 0x100000000ULL) >> 32);
        /* clang-format on */
        return ((x << 1) & 0xFFFFFFFDFFFFFFFFULL) | m;
      };

  auto multi_impl = [] [[using gnu: always_inline, const]] (
                        KMer::value_type x,
                        KMer::value_type n_rotations) noexcept {
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

namespace detail {

inline constexpr std::size_t kMaxK = 32;

inline constexpr auto kPrecomputed = [] consteval {
  std::array<std::array<std::uint64_t, kNtHashSeeds.size()>, kMaxK> dst;
  for (std::size_t i = 0; i < kMaxK; ++i) {
    for (std::size_t j = 0; j < kNtHashSeeds.size(); ++j) {
      dst[i][j] = srol(kNtHashSeeds[j], i);
    }
  }

  return dst;
}();

}  // namespace detail

inline constexpr auto nthash = [] [[using gnu: always_inline, pure, hot]]
                               (KMer::value_type prev, std::uint8_t base_out,
                                std::uint8_t base_in, std::uint64_t k) {
                                 return srol(prev) ^
                                        detail::kPrecomputed[k][base_out] ^
                                        kNtHashSeeds[base_in];
                               };
template <std::size_t N>
class alignas(64) Reg {
  std::int64_t data_[N];

 public:
  template <class Self>
  decltype(auto) operator[](this Self&& self, std::size_t i) {
    return std::forward_like<Self>(self.data_[i]);
  }
  std::int64_t* data() { return reinterpret_cast<std::int64_t*>(&data_); }
  std::size_t size() const noexcept { return 4; }
};

template <std::size_t N>
inline constexpr auto nthash_bulk = [] [[using gnu: always_inline, pure, hot]] (
                                        Reg<N> const& prev, Reg<N>& out,
                                        Reg<N>& in, std::uint64_t k) -> Reg<N> {
  Reg<N> dst;
  for (std::int64_t i = 0; i < prev.size(); ++i) {
    dst[i] = srol(prev[i]);
  }

  for (std::int64_t i = 0; i < out.size(); ++i) {
    out[i] = detail::kPrecomputed[k][out[i]];
  }

  for (std::int64_t i = 0; i < out.size(); ++i) {
    in[i] = kNtHashSeeds[in[i]];
  }

  *reinterpret_cast<__m256i*>(dst.data()) =
      _mm256_xor_si256(*reinterpret_cast<__m256i*>(dst.data()),
                       *reinterpret_cast<__m256i*>(out.data()));
  *reinterpret_cast<__m256i*>(dst.data()) =
      _mm256_xor_si256(*reinterpret_cast<__m256i*>(dst.data()),
                       *reinterpret_cast<__m256i*>(in.data()));

  return dst;
};

}  // namespace tb
