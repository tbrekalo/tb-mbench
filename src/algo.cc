#include "tb/algo.hpp"
#include "tb/ring.hpp"

#include <algorithm>
#include <cassert>
#include <deque>
#include <execution>
#include <span>

#include "eve/module/algo.hpp"

namespace tb {

namespace {

constexpr auto calc_mask =
    [] [[using gnu: always_inline, const]] (
        std::uint32_t kmer_length) constexpr noexcept -> std::uint64_t {
  return (1 << (kmer_length * 2)) - 1;
};

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

} // namespace

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

  auto push = [&front_idx, &back_idx,
               &dst](KMer::value_type hash_value,
                     KMer::position_type position) -> void {
    for (; front_idx < back_idx && dst[back_idx - 1].value() > hash_value;
         --back_idx)
      ;
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

std::vector<KMer> ArgminMinimize(MinimizeArgs args) {
  std::vector<KMer> dst;
  if (args.seq.size() < args.window_length + args.kmer_length - 2) {
    return dst;
  }

  dst.resize(args.seq.size());
  auto const mask = calc_mask(args.kmer_length);

  KMer::value_type value;
  std::vector<KMer::value_type> hashes;
  for (std::size_t i = 0; i < args.seq.size(); ++i) {
    value = ((value << 2) | args.seq.Code(i)) & mask;
    if (i >= args.kmer_length - 1) {
      hashes.push_back(hash(value, mask));
    }
  }

  std::size_t idx = 0, min_pos = 0;
  for (std::size_t i = args.window_length; i <= hashes.size(); ++i) {
    min_pos = std::min_element(hashes.begin() + i - args.window_length,
                               hashes.begin() + i) -
              hashes.begin();
    auto cond = idx == 0 || dst[idx - 1].position() != min_pos;
    min_pos = cond * min_pos + (!cond) * dst[idx].position();
    dst[idx] = KMer(hashes[min_pos], min_pos, 0);
    idx += cond;
  }

  dst.resize(idx);
  return dst;
}

std::vector<KMer> ArgminEveMinimize(MinimizeArgs args) {
  std::vector<KMer> dst;
  if (args.seq.size() < args.window_length + args.kmer_length - 2) {
    return dst;
  }

  dst.resize(args.seq.size());
  auto const mask = calc_mask(args.kmer_length);

  KMer::value_type value;
  std::vector<KMer::value_type> hashes;
  for (std::size_t i = 0; i < args.seq.size(); ++i) {
    value = ((value << 2) | args.seq.Code(i)) & mask;
    if (i >= args.kmer_length - 1) {
      hashes.push_back(hash(value, mask));
    }
  }

  std::int64_t idx = -1;
  for (std::size_t i = args.window_length; i <= hashes.size(); ++i) {
    auto window =
        std::span(hashes.begin() + i - args.window_length, hashes.begin() + i);
    if (auto min_pos = eve::algo::min_element(window) - window.begin() + i -
                       args.window_length;
        idx == -1 || dst[idx].position() != min_pos) {
      dst[++idx] = KMer(hashes[min_pos], min_pos, 0);
    }
  }

  dst.resize(idx + 1);
  return dst;
}

} // namespace tb
