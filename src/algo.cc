#include "tb/algo.hpp"

#include <deque>

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

struct WindowEntry {
  KMer::value_type hash_value;
  KMer::position_type position;
};

} // namespace

std::vector<KMer> NaiveMinimize(MinimizeArgs args) {
  std::vector<KMer> dst;
  if (args.seq.size() == 0) {
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
  if (args.seq.size() == 0) {
    return dst;
  }

  dst.reserve(args.seq.size());
  auto const mask = calc_mask(args.kmer_length);
  std::deque<WindowEntry> window;

  auto push = [&window](KMer::value_type hash_value,
                        KMer::position_type position) -> void {
    while (!window.empty() && window.back().hash_value > hash_value) {
      window.pop_back();
    }
    window.push_back(
        WindowEntry{.hash_value = hash_value, .position = position});
  };

  auto pop = [&window, w = args.window_length](KMer::position_type position) {
    if (window.front().position <= position - w) {
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
          (dst.empty() || dst.back().position() != window.front().position)) {
        dst.emplace_back(window.front().hash_value, window.front().position, 0);
      }
    }
  }

  return dst;
}

std::vector<KMer> InplaceMinimize(MinimizeArgs args) {
  std::vector<KMer> dst;
  if (args.seq.size() == 0) {
    return dst;
  }

  dst.reserve(args.seq.size());
  auto const mask = calc_mask(args.kmer_length);

  std::size_t front_idx = 0;
  std::size_t back_idx = 0;
  std::vector<KMer> window_buffer(args.seq.size());

  auto push = [&front_idx, &back_idx,
               &window_buffer](KMer::value_type hash_value,
                               KMer::position_type position) -> void {
    for (; front_idx < back_idx &&
           window_buffer[back_idx - 1].value() > hash_value;
         --back_idx)
      ;
    window_buffer[back_idx++] = KMer(hash_value, position, 0);
  };

  auto pop = [&front_idx, &window_buffer,
              w = args.window_length](KMer::position_type position) {
    if (window_buffer[front_idx].position() <= position - w) {
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
          (dst.empty() ||
           dst.back().position() != window_buffer[front_idx].position())) {
        dst.push_back(window_buffer[front_idx]);
      }
    }
  }

  return dst;
}

} // namespace tb
