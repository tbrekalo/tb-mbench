#include "tb/algo.hpp"

#include <algorithm>
#include <cassert>
#include <deque>
#include <span>

#include "eve/module/algo.hpp"
#include "tb/constants.hpp"
#include "tb/ring.hpp"

namespace tb {

namespace {

constexpr auto calc_mask =
    [] [[using gnu: always_inline, const]] (
        std::uint32_t kmer_length) constexpr noexcept -> std::uint64_t {
  return (1 << (kmer_length * 2)) - 1;
};

// Thomas Wang integer hash function
constexpr auto hash =
    [] [[using gnu: const, hot]] (
        std::uint64_t key,
        std::uint64_t mask) constexpr noexcept -> std::uint64_t {
  key = ((~key) + (key << 21)) & mask;
  key = key ^ (key >> 24);
  key = ((key + (key << 3)) + (key << 8)) & mask;
  key = key ^ (key >> 14);
  key = ((key + (key << 2)) + (key << 4)) & mask;
  key = key ^ (key >> 28);
  key = (key + (key << 31)) & mask;
  return key;
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

}  // namespace

std::vector<KMer> NaiveMinimize(MinimizeArgs args) {
  std::vector<KMer> dst;
  if (args.seq.size() < args.window_length + args.kmer_length - 2) {
    return dst;
  }

  dst.reserve(args.seq.size());
  auto const mask = calc_mask(args.kmer_length);
  for (std::size_t i = 0; i + args.kmer_length < args.seq.size(); ++i) {
    KMer::value_type min_hash;
    KMer::position_type min_position = args.seq.size();
    for (std::size_t j = 0; j < args.window_length &&
                            i + j + args.kmer_length - 1 < args.seq.size();
         ++j) {
      KMer::value_type value;
      for (std::size_t k = 0; k < args.kmer_length; ++k) {
        value = (value << 2) | args.seq.Code(i + j + k);
      }

      auto hash_value = hash(value, mask);
      if (min_position == args.seq.size() || hash_value < min_hash) {
        min_hash = hash_value;
        min_position = i + j;
      }
    }
    if (dst.empty() || dst.back().position() != min_position) {
      dst.emplace_back(min_hash, min_position, 0);
    }
  }

  return dst;
}

std::vector<KMer> DequeMinimize(MinimizeArgs args) {
  std::vector<KMer> dst;
  if (args.seq.size() < args.window_length + args.kmer_length - 2) {
    return dst;
  }

  dst.reserve(args.seq.size());
  auto const mask = calc_mask(args.kmer_length);
  std::deque<KMer> window;

  auto push = [&window](KMer::value_type hash_value,
                        KMer::position_type position) -> void {
    while (!window.empty() && window.back().value() > hash_value) {
      window.pop_back();
    }
    window.emplace_back(hash_value, position, 0);
  };

  auto pop = [&window, w = args.window_length](KMer::position_type position) {
    if (window.front().position() <= position - w) {
      window.pop_front();
    }
  };

  KMer::value_type value;
  for (std::size_t i = 0; i < args.seq.size(); ++i) {
    if (i >= args.window_length + args.kmer_length - 1) {
      pop(i - (args.kmer_length - 1));
    }

    value = ((value << 2) | args.seq.Code(i)) & mask;
    if (i >= args.kmer_length - 1) {
      push(hash(value, mask), i - (args.kmer_length - 1));
      if (i > args.window_length + args.kmer_length - 2 &&
          (dst.empty() || dst.back().position() != window.front().position())) {
        dst.emplace_back(window.front().value(), window.front().position(), 0);
      }
    }
  }

  return dst;
}

std::vector<KMer> InplaceMinimize(MinimizeArgs args) {
  std::vector<KMer> dst;
  if (args.seq.size() < args.window_length + args.kmer_length - 2) {
    return dst;
  }

  std::int64_t idx = -1;
  dst.resize(args.seq.size() - args.kmer_length);
  auto const mask = calc_mask(args.kmer_length);

  std::int64_t front_idx = 1;
  std::int64_t back_idx = 1;

  auto push = [&front_idx, &back_idx, &dst](
                  KMer::value_type hash_value,
                  KMer::position_type position) -> void {
    for (; front_idx < back_idx && dst[back_idx - 1].value() > hash_value;
         --back_idx);
    dst[back_idx++] = KMer(hash_value, position, 0);
  };

  auto pop = [&front_idx, &dst,
              w = args.window_length](KMer::position_type position) {
    if (dst[front_idx].position() <= position - w) {
      ++front_idx;
    }
  };

  KMer::value_type value;
  for (std::size_t i = 0; i < args.seq.size(); ++i) {
    if (i >= args.window_length + args.kmer_length - 1) {
      pop(i - (args.kmer_length - 1));
    }

    value = ((value << 2) | args.seq.Code(i)) & mask;
    if (i >= args.kmer_length - 1) {
      push(hash(value, mask), i - (args.kmer_length - 1));
      if (i > args.window_length + args.kmer_length - 2 &&
          (idx < 0 || dst[idx].position() != dst[front_idx].position())) {
        if (idx + 1 == front_idx) {
          std::shift_right(dst.begin() + front_idx++, dst.begin() + ++back_idx,
                           1);
        }
        assert(idx + 1 < front_idx);
        dst[++idx] = dst[front_idx];
      }
    }
  }

  dst.resize(idx + 1);
  return dst;
}

std::vector<KMer> RingMinimize(MinimizeArgs args) {
  std::vector<KMer> dst;
  if (args.seq.size() < args.window_length + args.kmer_length - 2) {
    return dst;
  }

  dst.reserve(args.seq.size());
  auto const mask = calc_mask(args.kmer_length);
  Ring window(args.window_length);

  auto push = [&window](KMer::value_type hash_value,
                        KMer::position_type position) -> void {
    while (!window.empty() && window.back().value() > hash_value) {
      window.pop_back();
    }
    window.push(KMer(hash_value, position, 0));
  };

  auto pop = [&window, w = args.window_length](KMer::position_type position) {
    if (window.front().position() <= position - w) {
      window.pop_front();
    }
  };

  KMer::value_type value;
  for (std::size_t i = 0; i < args.seq.size(); ++i) {
    if (i >= args.window_length + args.kmer_length - 1) {
      pop(i - (args.kmer_length - 1));
    }

    value = ((value << 2) | args.seq.Code(i)) & mask;
    if (i >= args.kmer_length - 1) {
      push(hash(value, mask), i - (args.kmer_length - 1));
      if (i > args.window_length + args.kmer_length - 2 &&
          (dst.empty() || dst.back().position() != window.front().position())) {
        dst.emplace_back(window.front().value(), window.front().position(), 0);
      }
    }
  }

  return dst;
}

namespace {

template <class Hasher, class Sampler>
class ArgMinMixinBase {
  [[no_unique_address]] Hasher hasher_;
  [[no_unique_address]] Sampler sampler_;

 public:
  std::vector<KMer> operator()(MinimizeArgs args) {
    if (args.seq.size() < args.window_length + args.kmer_length - 2) {
      return {};
    }

    return sampler_(args, hasher_(args));
  }
};

struct ThomasWangHasher {
  std::vector<KMer::value_type> operator()(MinimizeArgs args) const {
    auto const mask = calc_mask(args.kmer_length);

    KMer::value_type value;
    std::vector<KMer::value_type> hashes;
    for (std::size_t i = 0; i < args.seq.size(); ++i) {
      value = ((value << 2) | args.seq.Code(i)) & mask;
      if (i >= args.kmer_length - 1) {
        hashes.push_back(hash(value, mask));
      }
    }

    return hashes;
  }
};

struct NthHasher {
  std::vector<KMer::value_type> operator()(MinimizeArgs args) const {
    if (args.seq.size() < args.window_length + args.kmer_length - 2) {
      return {};
    }

    KMer::value_type value =
        srol(kNtHashSeeds[args.seq.Code(0)], args.kmer_length - 1);
    for (std::size_t i = 1; i < args.kmer_length; ++i) {
      value ^= srol(kNtHashSeeds[args.seq.Code(i)], args.kmer_length - (i + 1));
    }

    std::vector<KMer::value_type> hashes{value};
    for (std::size_t i = args.kmer_length; i < args.seq.size(); ++i) {
      value = nthash(value, args.seq.Code(i - (args.kmer_length - 1)),
                     args.seq.Code(i), args.kmer_length);
      hashes.push_back(value);
    }

    return hashes;
  }
};

struct ArgMinSampler {
  std::vector<KMer> operator()(MinimizeArgs args,
                               std::vector<KMer::value_type> hashes) const {
    std::vector<KMer> dst(hashes.size());
    std::int64_t idx = -1;
    for (std::size_t i = args.window_length; i <= hashes.size(); ++i) {
      if (auto min_pos =
              std::min_element(hashes.begin() + i - args.window_length,
                               hashes.begin() + i) -
              hashes.begin();
          idx == -1 || dst[idx].position() != min_pos) {
        dst[++idx] = KMer(hashes[min_pos], min_pos, 0);
      }
    }

    dst.resize(idx + 1);
    return dst;
  }
};

struct ArgMinRecoverySampler {
  std::vector<KMer> operator()(MinimizeArgs args,
                               std::vector<KMer::value_type> hashes) {
    std::vector<KMer> dst(hashes.size());
    std::size_t min_pos =
        std::min_element(hashes.begin(), hashes.begin() + args.window_length) -
        hashes.begin();
    dst[0] = KMer(hashes[min_pos], min_pos, 0);

    std::size_t idx = 1;
    for (std::size_t i = args.window_length + 1; i <= hashes.size(); ++i) {
      bool cond;
      if (min_pos >= i - args.window_length) {
        cond = hashes[i - 1] < hashes[min_pos];
        min_pos = cond * (i - 1) + (!cond) * min_pos;
      } else {
        min_pos = std::min_element(hashes.begin() + i - args.window_length,
                                   hashes.begin() + i) -
                  hashes.begin();
        cond = dst[idx - 1].position() != min_pos;
        min_pos = cond * min_pos + (!cond) * dst[idx].position();
      }
      dst[idx] = KMer(hashes[min_pos], min_pos, 0);
      idx += cond;
    }

    dst.resize(idx);
    return dst;
  }
};

struct ArgMinEveRecoverySampler {
  std::vector<KMer> operator()(MinimizeArgs args,
                               std::vector<KMer::value_type> hashes) {
    std::vector<KMer> dst(hashes.size());
    std::size_t min_pos =
        std::min_element(hashes.begin(), hashes.begin() + args.window_length) -
        hashes.begin();
    dst[0] = KMer(hashes[min_pos], min_pos, 0);

    std::size_t idx = 1;
    for (std::size_t i = args.window_length + 1; i <= hashes.size(); ++i) {
      bool cond;
      if (min_pos >= i - args.window_length) {
        cond = hashes[i - 1] < hashes[min_pos];
        min_pos = cond * (i - 1) + (!cond) * min_pos;
      } else {
        auto window = std::span(hashes.begin() + i - args.window_length,
                                hashes.begin() + i);
        min_pos = eve::algo::min_element(window) - window.begin() + i -
                  args.window_length;
        cond = dst[idx - 1].position() != min_pos;
        min_pos = cond * min_pos + (!cond) * dst[idx].position();
      }
      dst[idx] = KMer(hashes[min_pos], min_pos, 0);
      idx += cond;
    }

    dst.resize(idx);
    return dst;
  }
};

using ArgMinMixin = ArgMinMixinBase<ThomasWangHasher, ArgMinSampler>;
using NtHashArgMinMixin = ArgMinMixinBase<NthHasher, ArgMinSampler>;
using ArgMinRecoveryMixin =
    ArgMinMixinBase<ThomasWangHasher, ArgMinRecoverySampler>;
using ArgMinEverRecoveryMixin =
    ArgMinMixinBase<ThomasWangHasher, ArgMinEveRecoverySampler>;
using NtHashArgMinRecoveryMixin =
    ArgMinMixinBase<NthHasher, ArgMinRecoverySampler>;

}  // namespace

std::vector<KMer> ArgMinMinimize(MinimizeArgs args) {
  return ArgMinMixin{}(args);
}

std::vector<KMer> NtHashArgMinMinimize(MinimizeArgs args) {
  return NtHashArgMinMixin{}(args);
}

std::vector<KMer> ArgMinRecoveryMinimize(MinimizeArgs args) {
  return ArgMinRecoveryMixin{}(args);
}

std::vector<KMer> ArgMinRecoveryEveMinimize(MinimizeArgs args) {
  return ArgMinEverRecoveryMixin{}(args);
}

std::vector<KMer> NtHashArgMinRecoveryMinimize(MinimizeArgs args) {
  return NtHashArgMinRecoveryMixin{}(args);
}

}  // namespace tb
