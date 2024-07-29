#include "tb/algo.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <concepts>
#include <deque>
#include <span>

#include "eve/module/algo.hpp"
#include "tb/nthash.hpp"
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
  std::vector<KMer> operator()(MinimizeArgs args) const {
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
    std::vector<KMer::value_type> hashes(args.seq.size() - args.kmer_length +
                                         1);
    for (std::size_t i = 0; i < args.seq.size(); ++i) {
      value = ((value << 2) | args.seq.Code(i)) & mask;
      if (i >= args.kmer_length - 1) {
        hashes[i - (args.kmer_length - 1)] = hash(value, mask);
      }
    }

    return hashes;
  }
};

template <auto roll_fn>
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

    std::vector<KMer::value_type> hashes(args.seq.size() - args.kmer_length +
                                         1);
    hashes[0] = value;

    for (std::size_t i = args.kmer_length; i < args.seq.size(); ++i) {
      value = roll_fn(value, args.seq.Code(i - args.kmer_length),
                      args.seq.Code(i), args.kmer_length);
      hashes[i - args.kmer_length] = value;
    }

    return hashes;
  }
};

template <class T>
concept AMinElement = requires(T lhs, T rhs) {
  { lhs < rhs } -> std::same_as<bool>;
  requires std::copyable<T>;
};

struct StdMinElement {
  template <AMinElement T>
  constexpr std::span<T>::iterator operator()(
      std::span<T> span) const noexcept {
    return std::ranges::min_element(span);
  }
};

struct EveMinElement {
  template <AMinElement T>
  constexpr std::span<T>::iterator operator()(
      std::span<T> span) const noexcept {
    return eve::algo::min_element(span);
  }
};

struct PredicationMinElement {
  template <AMinElement T>
  constexpr std::span<T>::iterator operator()(
      std::span<T> span) const noexcept {
    auto idx = 0;
    for (std::size_t jdx = 0; jdx < span.size(); ++jdx) {
      auto c = span[jdx] < span[idx];
      idx = c * jdx + (1 - c) * idx;
    }

    return span.begin() + idx;
  }
};

struct DuffMinElement {
  template <AMinElement T>
  constexpr std::span<T>::iterator operator()(
      std::span<T> span) const noexcept {
    std::int64_t i = span.size(), min_pos = 0;
    auto min = [span](std::int64_t lhs_idx, std::int64_t rhs_idx) {
      auto cond = span[lhs_idx] <= span[rhs_idx];
      return cond * lhs_idx + (1 - cond) * rhs_idx;
    };

    switch (i % 4) {
      while (i > 0) {
        case 0:
          min_pos = min(--i, min_pos);
        case 3:
          min_pos = min(--i, min_pos);
        case 2:
          min_pos = min(--i, min_pos);
        case 1:
          min_pos = min(--i, min_pos);
      }
    }

    return span.begin() + min_pos;
  }
};

template <class MinPolicy>
class ArgMinSampler {
  [[no_unique_address]] MinPolicy min_element_;

 public:
  std::vector<KMer> operator()(
      MinimizeArgs args, std::vector<KMer::value_type> hashes) const noexcept {
    std::vector<KMer> dst(hashes.size());
    std::int64_t idx = -1;
    for (std::size_t i = args.window_length; i <= hashes.size(); ++i) {
      auto window = std::span(hashes.begin() + i - args.window_length,
                              hashes.begin() + i);
      if (auto min_pos =
              min_element_(window) - window.begin() + i - args.window_length;
          idx == -1 || dst[idx].position() != min_pos) {
        dst[++idx] = KMer(hashes[min_pos], min_pos, 0);
      }
    }

    dst.resize(idx + 1);
    return dst;
  }
};

template <class MinPolicy>
class ArgMinRecoverySampler {
  [[no_unique_address]] MinPolicy min_element_;

 public:
  std::vector<KMer> operator()(
      MinimizeArgs args, std::vector<KMer::value_type> hashes) const noexcept {
    std::vector<KMer> dst(hashes.size());
    std::size_t min_pos =
        std::min_element(hashes.begin(), hashes.begin() + args.window_length) -
        hashes.begin();
    dst[0] = KMer(hashes[min_pos], min_pos, 0);

    std::size_t idx = 1;
    for (std::size_t i = args.window_length + 1; i <= hashes.size(); ++i) {
      bool cond = 1;
      if (min_pos >= i - args.window_length) {
        cond = hashes[i - 1] < hashes[min_pos];
        min_pos = cond * (i - 1) + (!cond) * min_pos;
      } else {
        auto window = std::span(hashes.begin() + i - args.window_length,
                                hashes.begin() + i);
        min_pos =
            min_element_(window) - window.begin() + i - args.window_length;
      }
      dst[idx] = KMer(hashes[min_pos], min_pos, 0);
      idx += cond;
    }

    dst.resize(idx);
    return dst;
  }
};

template <template <class> class Sampler>
class UnrolledSampler {
  static constexpr std::size_t kMaxW = 31;
  static constexpr std::size_t kJumpTblSize = kMaxW + 2uz;

  using ImplPtr = std::vector<KMer> (*)(MinimizeArgs,
                                        std::vector<KMer::value_type>);

  template <std::size_t I>
  static constexpr auto ImplGenerator = []() -> ImplPtr {
    return +[](MinimizeArgs args,
               std::vector<KMer::value_type> hashes) -> std::vector<KMer> {
      return Sampler<decltype([] [[using gnu: always_inline, hot, const]] (
                                  std::span<KMer::value_type> span) {
        auto min = 0;
        [&]<std::size_t... Is>(std::index_sequence<Is...>) {
          (..., [&](auto j) {
            auto c = span[j] < span[min];
            min = c * j + (1 - c) * min;
          }(Is));
        }(std::make_index_sequence<I>{});
        return span.begin() + min;
      })>{}(args, std::move(hashes));
    };
  };

  static constexpr auto kJumpTable = [] consteval {
    std::array<ImplPtr, kJumpTblSize> dst;
    [&]<std::size_t... Is>(std::index_sequence<Is...>) {
      (..., (dst[Is] = ImplGenerator<Is>()));
    }(std::make_index_sequence<kJumpTblSize>{});

    return dst;
  }();

 public:
  std::vector<KMer> operator()(
      MinimizeArgs args, std::vector<KMer::value_type> hashes) const noexcept {
    return kJumpTable[args.window_length](args, std::move(hashes));
  }
};

struct ArgMinRolling {
  std::vector<KMer> operator()(MinimizeArgs args,
                               std::vector<KMer::value_type> hashes) const {
    std::vector<KMer> dst(hashes.size());
    std::int64_t idx = 0, min_pos = std::min_element(
                                        hashes.begin(),
                                        hashes.begin() + args.window_length) -
                                    hashes.begin();

    dst[idx++] = KMer(hashes[min_pos], min_pos, 0);
    min_pos = std::min_element(hashes.begin() + min_pos + 1,
                               hashes.begin() + args.window_length) -
              hashes.begin();

    for (std::int64_t i = args.window_length; i < hashes.size(); ++i) {
      min_pos = hashes[i] < hashes[min_pos] ? i : min_pos;
      if (dst[idx - 1].position() > i - args.window_length &&
          hashes[min_pos] >= dst[idx - 1].value()) {
        continue;
      }
      dst[idx++] = KMer(hashes[min_pos], min_pos, 0);
      i = min_pos++;
    }

    dst.resize(idx);
    return dst;
  }
};

// Initialize ArgMin samplers
using StdArgMinSampler = ArgMinSampler<PredicationMinElement>;
using EveArgMinSampler = ArgMinSampler<EveMinElement>;
using PredicationArgMinSampler = ArgMinSampler<PredicationMinElement>;
using DuffArgMinSampler = ArgMinSampler<DuffMinElement>;
using UnrolledArgMinSampler = UnrolledSampler<ArgMinSampler>;

// Initialize ArgMinRecovery samplers
using StdArgMinRecoverySampler = ArgMinRecoverySampler<PredicationMinElement>;
using EveArgMinRecoverySampler = ArgMinRecoverySampler<EveMinElement>;
using DuffArgMinRecoverySampler = ArgMinRecoverySampler<DuffMinElement>;
using PredicationArgMinRecoverySampler =
    ArgMinRecoverySampler<PredicationMinElement>;
using UnrolledArgMinRecoverySampler = UnrolledSampler<ArgMinRecoverySampler>;

// ArgMin mixins
using ArgMinMixin = ArgMinMixinBase<ThomasWangHasher, PredicationArgMinSampler>;
using ArgMinEveMixin = ArgMinMixinBase<ThomasWangHasher, EveArgMinSampler>;
using ArgMinDuffMixin = ArgMinMixinBase<ThomasWangHasher, DuffArgMinSampler>;
using ArgMinUnrolledMixin =
    ArgMinMixinBase<ThomasWangHasher, UnrolledArgMinSampler>;
// NtHash ArgMin mixins
using NtHashArgMinMixin =
    ArgMinMixinBase<NthHasher<nthash<NtHashImpl::kRuntime>>,
                    PredicationArgMinSampler>;
using NtHashPrecomputedArgMinMixin =
    ArgMinMixinBase<NthHasher<nthash<NtHashImpl::kPrecomputed>>,
                    PredicationArgMinSampler>;
using NtHashPrecomputedArgMinUnrolledMixin =
    ArgMinMixinBase<NthHasher<nthash<NtHashImpl::kPrecomputed>>,
                    UnrolledArgMinSampler>;

// ArgMinRecovery mixins
using ArgMinRecoveryMixin =
    ArgMinMixinBase<ThomasWangHasher, PredicationArgMinRecoverySampler>;
using ArgMinEveRecoveryMixin =
    ArgMinMixinBase<ThomasWangHasher, EveArgMinRecoverySampler>;
using ArgMinDuffRecoveryMixin =
    ArgMinMixinBase<ThomasWangHasher, DuffArgMinRecoverySampler>;
using ArgMinUnrolledRecoveryMixin =
    ArgMinMixinBase<ThomasWangHasher, UnrolledArgMinRecoverySampler>;
// NtHash ArgMin recovery mixins
using NtHashArgMinRecoveryMixin =
    ArgMinMixinBase<NthHasher<nthash<NtHashImpl::kRuntime>>,
                    PredicationArgMinRecoverySampler>;
using NtHashPrecomputedArgMinRecoveryMixin =
    ArgMinMixinBase<NthHasher<nthash<NtHashImpl::kPrecomputed>>,
                    PredicationArgMinRecoverySampler>;
using NtHashPrecomputedArgMinUnrolledRecoveryMixin =
    ArgMinMixinBase<NthHasher<nthash<NtHashImpl::kPrecomputed>>,
                    UnrolledArgMinRecoverySampler>;

// ArgMinRolling mixins
using ArgMinRollingMixin = ArgMinMixinBase<ThomasWangHasher, ArgMinRolling>;
}  // namespace

// Arg min based implementations
std::vector<KMer> ArgMinMinimize(MinimizeArgs args) {
  return ArgMinMixin{}(args);
}

std::vector<KMer> ArgMinEveMinimize(MinimizeArgs args) {
  return ArgMinEveMixin{}(args);
}

std::vector<KMer> ArgMinDuffMinimize(MinimizeArgs args) {
  return ArgMinDuffMixin{}(args);
}

std::vector<KMer> ArgMinUnrolledMinimize(MinimizeArgs args) {
  return ArgMinUnrolledMixin{}(args);
}

std::vector<KMer> NtHashArgMinMinimize(MinimizeArgs args) {
  return NtHashArgMinMixin{}(args);
}

std::vector<KMer> NtHashPrecomputedArgMinMinimize(MinimizeArgs args) {
  return NtHashPrecomputedArgMinMixin{}(args);
}

std::vector<KMer> NtHashPrecomputedArgMinUnrolledMinimize(MinimizeArgs args) {
  return NtHashPrecomputedArgMinUnrolledMixin{}(args);
}

// Arg min recovery based implementations
std::vector<KMer> ArgMinRecoveryMinimize(MinimizeArgs args) {
  return ArgMinRecoveryMixin{}(args);
}

std::vector<KMer> ArgMinRecoveryEveMinimize(MinimizeArgs args) {
  return ArgMinEveRecoveryMixin{}(args);
}

std::vector<KMer> ArgMinRecoveryDuffMinimize(MinimizeArgs args) {
  return ArgMinDuffRecoveryMixin{}(args);
}

std::vector<KMer> ArgMinRecoveryUnrolledMinimize(MinimizeArgs args) {
  return ArgMinUnrolledRecoveryMixin{}(args);
}

std::vector<KMer> NtHashArgMinRecoveryMinimize(MinimizeArgs args) {
  return NtHashArgMinRecoveryMixin{}(args);
}

std::vector<KMer> NtHashPrecomputedArgMinRecoveryMinimize(MinimizeArgs args) {
  return NtHashPrecomputedArgMinRecoveryMixin{}(args);
}

std::vector<KMer> NtHashPrecomputedArgMinUnrolledRecoveryMinimize(
    MinimizeArgs args) {
  return NtHashPrecomputedArgMinUnrolledRecoveryMixin{}(args);
}

std::vector<KMer> ArgMinRollingMinimize(MinimizeArgs args) {
  return ArgMinRollingMixin{}(args);
}

}  // namespace tb
