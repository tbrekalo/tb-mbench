#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "gtest/gtest.h"

#include "tb/algo.hpp"

namespace {

constexpr int kSeed = 42;

class MinimizeTest : public testing::Test {
protected:
  MinimizeTest()
      : seq_(1uz << 14uz, kSeed), args_(tb::MinimizeArgs{
                                      .seq = seq_,
                                      .window_length = 5,
                                      .kmer_length = 15,
                                  }) {}

  tb::MockSequence seq_;
  tb::MinimizeArgs args_;
};

} // namespace

TEST_F(MinimizeTest, NaiveDensity) {
  auto minimizers = tb::NaiveMinimize(args_);
  EXPECT_GE(minimizers.size(),
            static_cast<double>(seq_.size()) / args_.window_length);
}

TEST_F(MinimizeTest, DequeVsNaive) {
  auto naive_minimizers = tb::NaiveMinimize(args_);
  auto deque_minimizers = tb::DequeMinimize(args_);

  EXPECT_EQ(naive_minimizers, deque_minimizers);
}

TEST_F(MinimizeTest, InplaceVsNaive) {
  auto naive_minimizers = tb::NaiveMinimize(args_);
  auto inplace_minimizers = tb::InplaceMinimize(args_);

  EXPECT_EQ(naive_minimizers, inplace_minimizers);
}

TEST_F(MinimizeTest, RingVsNaive) {
  auto naive_minimizers = tb::NaiveMinimize(args_);
  auto ring_minimizers = tb::RingMinimize(args_);

  EXPECT_EQ(naive_minimizers, ring_minimizers);
}
