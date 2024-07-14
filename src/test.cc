#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "gtest/gtest.h"

#include "tb/algo.hpp"

namespace {

constexpr int kSeed = 42;

class MinimizeTest : public testing::Test {
protected:
  MinimizeTest()
      : seq_(1uz << 14uz, kSeed), kmer_length_(15), window_length_(5) {}

  tb::MockSequence seq_;
  std::uint32_t kmer_length_;
  std::uint32_t window_length_;
};

} // namespace

TEST_F(MinimizeTest, Density) {
  auto minimizers = tb::NaiveMinimize({
      .seq = seq_,
      .kmer_length = kmer_length_,
      .window_length = window_length_,
  });
  EXPECT_TRUE(minimizers.size() >=
              static_cast<double>(seq_.n_bases()) / window_length_);
}
