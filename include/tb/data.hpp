#pragma once

#include <cstdint>
#include <vector>

namespace tb {

// biosoup::NucleicAcid like structure
class MockSequence {
  std::size_t n_bases_;
  std::vector<std::uint64_t> data_;

public:
  MockSequence(std::size_t n_bases, int seed);

  [[gnu::always_inline]] std::uint64_t Code(std::size_t i) const noexcept {
    return ((data_[i >> 5] >> ((i << 1) & 63)) & 3);
  }

  [[gnu::always_inline]] std::uint64_t
  ReverseCode(std::size_t i) const noexcept {
    return Code(n_bases_ - i - 1) ^ 3;
  }

  std::size_t size() const noexcept { return n_bases_; }
};

class KMer {
  using ValueTypeImpl = std::uint64_t;

  ValueTypeImpl value_;
  std::uint64_t pos_strand_;

public:
  using value_type = ValueTypeImpl;
  using position_type = std::int32_t;

  KMer() = default;
  KMer(value_type value, position_type pos, bool strand)
      : value_(value),
        pos_strand_((static_cast<std::uint64_t>(strand) << 31) | pos) {}

  [[gnu::always_inline]] value_type value() const noexcept { return value_; }
  [[gnu::always_inline]] position_type position() const noexcept {
    return pos_strand_ & ((1ul << 31ul) - 1ul);
  }
  [[gnu::always_inline]] bool strand() const noexcept {
    return (pos_strand_ >> 31) & 1;
  }

  friend auto operator<=>(KMer const &lhs, KMer const &rhs) = default;
};

} // namespace tb
