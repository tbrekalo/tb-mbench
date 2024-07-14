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
            static_cast<double>(seq_.n_bases()) / args_.window_length);
}

TEST_F(MinimizeTest, DequeAgainsNaive) {
  auto naive_minimizers = tb::NaiveMinimize(args_);
  auto deque_minimizers = tb::DequeMinimize(args_);

  EXPECT_EQ(naive_minimizers, deque_minimizers);
}
