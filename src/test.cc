#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "gtest/gtest.h"
#include "tb/algo.hpp"

namespace {

constexpr int kSeed = 42;

class MinimizeTest : public testing::Test {
 protected:
  MinimizeTest()
      : seq_(1uz << 14uz, kSeed),
        args_(tb::MinimizeArgs{
            .seq = seq_,
            .window_length = 5,
            .kmer_length = 15,
        }) {}

  tb::MockSequence seq_;
  tb::MinimizeArgs args_;
};

}  // namespace

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

TEST_F(MinimizeTest, ArgMinVsNaive) {
  auto naive_minimizers = tb::NaiveMinimize(args_);
  auto argmin_minimizers = tb::ArgMinMinimize(args_);

  EXPECT_EQ(naive_minimizers, argmin_minimizers);
}

TEST_F(MinimizeTest, ArgMinUnrolledVsArgMin) {
  auto argmin_minimizers = tb::ArgMinMinimize(args_);
  auto unrolled_minimizers = tb::ArgMinUnrolledMinimize(args_);

  EXPECT_EQ(argmin_minimizers, unrolled_minimizers);
}

TEST_F(MinimizeTest, ArgMinRecoveryVsNaive) {
  auto naive_minimizers = tb::NaiveMinimize(args_);
  auto argmin_minimizers = tb::ArgMinRecoveryMinimize(args_);

  EXPECT_EQ(naive_minimizers, argmin_minimizers);
}

TEST_F(MinimizeTest, ArgMinRecoveryUnrolledVsNaive) {
  auto naive_minimizers = tb::NaiveMinimize(args_);
  auto argmin_minimizers = tb::ArgMinRecoveryUnrolledMinimize(args_);

  EXPECT_EQ(naive_minimizers, argmin_minimizers);
}


TEST_F(MinimizeTest, SplitWindowVsArgMin) {
  auto argmin_minimizers = tb::ArgMinMinimize(args_);
  auto split_minimizers = tb::SplitWindowMinimize(args_);

  EXPECT_EQ(argmin_minimizers, split_minimizers);
}
