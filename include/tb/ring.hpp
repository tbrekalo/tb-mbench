#pragma once

#include "tb/data.hpp"

namespace tb {

class Ring {
  std::int64_t front_idx_;
  std::int64_t back_idx_;
  std::vector<KMer> data_;

public:
  Ring(std::size_t size) : front_idx_(0), back_idx_(0), data_(size) {}

  KMer const &front() const noexcept { return data_[front_idx_]; }

  KMer const &back() const noexcept {
    return data_[(data_.size() + back_idx_ - 1) % data_.size()];
  }

  bool empty() const noexcept { return front_idx_ == back_idx_; }

  void push(KMer kmer) {
    data_[back_idx_] = kmer;
    back_idx_ = (back_idx_ + 1) % data_.size();
  }

  void pop_back() noexcept {
    back_idx_ = (data_.size() + back_idx_ - 1) % data_.size();
  }

  void pop_front() noexcept { front_idx_ = (front_idx_ + 1) % data_.size(); }
};

} // namespace tb
