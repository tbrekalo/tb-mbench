#include "tb/algo.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <concepts>
#include <deque>
#include <span>

#include "tb/nthash.hpp"

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

struct NtHasher {
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

    for (std::int64_t i = args.kmer_length; i < args.seq.size(); ++i) {
      value = nthash(value, args.seq.Code(i - args.kmer_length),
                     args.seq.Code(i), args.kmer_length);
      hashes[i - args.kmer_length + 1] = value;
    }

    return hashes;
  }
};

class NtHasherOpt {
  template <std::size_t N>
  std::vector<KMer::value_type> impl(MinimizeArgs args) const {
    if (args.seq.size() < N * args.kmer_length) {
      return NtHasher{}(args);
    }

    std::int64_t n_kmers = args.seq.size() - args.kmer_length + 1;
    std::vector<KMer::value_type> dst(n_kmers);

    auto pivots = [n_kmers] {
      std::array<std::int64_t, N> pivots;
      std::int64_t pivot = n_kmers / N;

      for (std::int64_t i = 0; i < N; ++i) {
        pivots[i] = pivot * (i + 1);
      }

      return pivots;
    }();

    auto indices = [n_kmers] {
      std::array<std::int64_t, N> indices;
      std::int64_t pivot = n_kmers / N;

      for (std::int64_t i = 0; i < N; ++i) {
        indices[i] = pivot * i;
      }

      return indices;
    }();

    auto idx = indices;
    auto values = [args, &indices] {
      std::array<KMer::value_type, N> values;
      for (std::int64_t i = 0; i < N; ++i) {
        values[i] = 0;
      }

      for (std::int64_t i = 0; i < args.kmer_length; ++i) {
        for (std::int64_t j = 0; j < N; ++j) {
          values[j] ^= srol(kNtHashSeeds[args.seq.Code(indices[j]++)],
                            args.kmer_length - (i + 1));
        }
      }
      return values;
    }();

    for (std::int64_t i = 0; i < N; ++i) {
      dst[idx[i]++] = values[i];
    }

    using RegType = Reg<N>;
    while (indices.front() - args.kmer_length + 1 < pivots.front()) {
      RegType v, base_out, base_in;
      for (std::int64_t i = 0; i < N; ++i) {
        v[i] = values[i];
        base_out[i] = args.seq.Code(indices[i] - args.kmer_length);
        base_in[i] = args.seq.Code(indices[i]);
      }

      v = nthash_bulk<N>(v, base_out, base_in, args.kmer_length);
      for (std::int64_t i = 0; i < N; ++i) {
        values[i] = v[i];
        dst[idx[i]++] = values[i];
        ++indices[i];
      }
    }

    for (; indices.back() < args.seq.size(); ++indices.back()) {
      values.back() = nthash(values.back(),
                             args.seq.Code(indices.back() - args.kmer_length),
                             args.seq.Code(indices.back()), args.kmer_length);
      dst[idx.back()++] = values.back();
    }

    return dst;
  }

 public:
  std::vector<KMer::value_type> operator()(MinimizeArgs args) const {
    return impl<4>(args);
  }
};

template <class T>
concept AMinElement = requires(T lhs, T rhs) {
  { lhs < rhs } -> std::same_as<bool>;
  requires std::copyable<T>;
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

class SplitWindow {
  std::vector<KMer> impl(MinimizeArgs args,
                         std::vector<KMer::value_type> hashes) const {
    if (hashes.empty()) {
      return {};
    }

    std::vector<KMer> dst(hashes.size());
    std::int64_t idx = 0;

    std::vector<std::int64_t> lhs(args.window_length + 1);
    std::int64_t lhs_idx = 1;

    std::vector<std::int64_t> rhs(args.window_length + 1);
    std::int64_t rhs_idx = 1, rhs_min = 0;

    auto shift_stacks = [&] {
      for (; rhs_idx > 1; --rhs_idx) {
        auto cond = lhs_idx == 1 ||
                    hashes[rhs[rhs_idx - 1]] <= hashes[lhs[lhs_idx - 1]];
        lhs[lhs_idx++] =
            cond * rhs[rhs_idx - 1] + (1 - cond) * lhs[lhs_idx - 1];
      }
    };

    auto push_back = [&](std::int64_t i) {
      rhs_min = hashes[i] < hashes[rhs_min] ? i : rhs_min;
      rhs[rhs_idx++] = i;
    };

    auto pop_front = [&](std::int64_t i) {
      if (lhs_idx == 1) {
        shift_stacks();
        rhs_min = i;
      }
      --lhs_idx;
    };

    for (std::int64_t i = 0;
         i < std::min<std::int64_t>(hashes.size(), args.window_length); ++i) {
      push_back(i);
    }

    dst[idx++] = KMer(hashes[rhs_min], rhs_min, 0);
    pop_front(args.window_length);
    for (std::int64_t i = args.window_length; i < hashes.size();
         i += args.window_length) {
      for (std::int64_t j = 0; j < args.window_length && i + j < hashes.size();
           ++j) {
        push_back(i + j);
        auto min_pos =
            lhs_idx > 1 && hashes[lhs[lhs_idx - 1]] <= hashes[rhs_min]
                ? lhs[lhs_idx - 1]
                : rhs_min;
        if (dst[idx - 1].position() != min_pos) {
          dst[idx++] = KMer(hashes[min_pos], min_pos, 0);
        }

        pop_front(i + j + 1);
      }
    }

    dst.resize(idx);
    return dst;
  }

 public:
  std::vector<KMer> operator()(MinimizeArgs args,
                               std::vector<KMer::value_type> hashes) const {
    return impl(args, hashes);
  }
};

// Initialize ArgMin samplers
using PredicationArgMinSampler = ArgMinSampler<PredicationMinElement>;
using UnrolledArgMinSampler = UnrolledSampler<ArgMinSampler>;

// Initialize ArgMinRecovery samplers
using PredicationArgMinRecoverySampler =
    ArgMinRecoverySampler<PredicationMinElement>;
using UnrolledArgMinRecoverySampler = UnrolledSampler<ArgMinRecoverySampler>;

// ArgMin mixins
using ArgMinMixin = ArgMinMixinBase<ThomasWangHasher, PredicationArgMinSampler>;
using ArgMinUnrolledMixin =
    ArgMinMixinBase<ThomasWangHasher, UnrolledArgMinSampler>;
// NtHash ArgMin mixins
using NtHashPrecomputedArgMinUnrolledMixin =
    ArgMinMixinBase<NtHasher, UnrolledArgMinSampler>;

// ArgMinRecovery mixins
using ArgMinRecoveryMixin =
    ArgMinMixinBase<ThomasWangHasher, PredicationArgMinRecoverySampler>;
using ArgMinUnrolledRecoveryMixin =
    ArgMinMixinBase<ThomasWangHasher, UnrolledArgMinRecoverySampler>;
// NtHash ArgMin recovery mixins
using NtHashPrecomputedArgMinUnrolledRecoveryMixin =
    ArgMinMixinBase<NtHasher, UnrolledArgMinRecoverySampler>;

// SplitWindow mixins
using SplitWindowMixin = ArgMinMixinBase<ThomasWangHasher, SplitWindow>;

}  // namespace

std::vector<KMer::value_type> NtHash(MinimizeArgs args) {
  return NtHasher{}(args);
}

std::vector<KMer::value_type> NtHashOpt(MinimizeArgs args) {
  return NtHasherOpt{}(args);
}

// Arg min based implementations
std::vector<KMer> ArgMinMinimize(MinimizeArgs args) {
  return ArgMinMixin{}(args);
}

std::vector<KMer> ArgMinUnrolledMinimize(MinimizeArgs args) {
  return ArgMinUnrolledMixin{}(args);
}

std::vector<KMer> NtHashArgMinMinimize(MinimizeArgs args) {
  return NtHashPrecomputedArgMinUnrolledMixin{}(args);
}

// Arg min recovery based implementations
std::vector<KMer> ArgMinRecoveryMinimize(MinimizeArgs args) {
  return ArgMinRecoveryMixin{}(args);
}

std::vector<KMer> ArgMinRecoveryUnrolledMinimize(MinimizeArgs args) {
  return ArgMinUnrolledRecoveryMixin{}(args);
}

std::vector<KMer> NtHashRecoveryMinimize(MinimizeArgs args) {
  return NtHashPrecomputedArgMinUnrolledRecoveryMixin{}(args);
}

std::vector<KMer> SplitWindowMinimize(MinimizeArgs args) {
  return SplitWindowMixin{}(args);
}

}  // namespace tb
