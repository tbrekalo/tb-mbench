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

TEST_F(MinimizeTest, ArgMinVsNaive) {
  auto naive_minimizers = tb::NaiveMinimize(args_);
  auto argmin_minimizers = tb::ArgMinMinimize(args_);

  EXPECT_EQ(naive_minimizers, argmin_minimizers);
}

TEST_F(MinimizeTest, ArgMinEveVsNaive) {
  auto naive_minimizers = tb::NaiveMinimize(args_);
  auto argmin_minimizers = tb::ArgMinEveMinimize(args_);

  EXPECT_EQ(naive_minimizers, argmin_minimizers);
}

TEST_F(MinimizeTest, ArgMinDuffVsNaive) {
  auto naive_minimizers = tb::NaiveMinimize(args_);
  auto argmin_minimizers = tb::ArgMinDuffMinimize(args_);

  EXPECT_EQ(naive_minimizers, argmin_minimizers);
}

TEST_F(MinimizeTest, ArgMinUnrolledVsNaive) {
  auto naive_minimizers = tb::NaiveMinimize(args_);
  auto argmin_minimizers = tb::ArgMinUnrolledMinimize(args_);

  EXPECT_EQ(naive_minimizers, argmin_minimizers);
}

TEST_F(MinimizeTest, ArgMinRecoveryVsNaive) {
  auto naive_minimizers = tb::NaiveMinimize(args_);
  auto argmin_minimizers = tb::ArgMinRecoveryMinimize(args_);

  EXPECT_EQ(naive_minimizers, argmin_minimizers);
}

TEST_F(MinimizeTest, ArgMinRecoveryEveVsNaive) {
  auto naive_minimizers = tb::NaiveMinimize(args_);
  auto argmin_minimizers = tb::ArgMinRecoveryEveMinimize(args_);

  EXPECT_EQ(naive_minimizers, argmin_minimizers);
}

TEST_F(MinimizeTest, ArgMinRecoveryDuffVsNaive) {
  auto naive_minimizers = tb::NaiveMinimize(args_);
  auto argmin_minimizers = tb::ArgMinRecoveryDuffMinimize(args_);

  EXPECT_EQ(naive_minimizers, argmin_minimizers);
}

TEST_F(MinimizeTest, ArgMinRecoveryUnrolledVsNaive) {
  auto naive_minimizers = tb::NaiveMinimize(args_);
  auto argmin_minimizers = tb::ArgMinRecoveryUnrolledMinimize(args_);

  EXPECT_EQ(naive_minimizers, argmin_minimizers);
}

TEST_F(MinimizeTest, NtHashPrecomputedArgMinVsNtHashArgMin) {
  auto minimizers = tb::NtHashArgMinMinimize(args_);
  auto precomputed_minimizers = tb::NtHashPrecomputedArgMinMinimize(args_);

  EXPECT_EQ(minimizers, precomputed_minimizers);
}

TEST_F(MinimizeTest, NtHashPrecomputedArgMinUnrolledVsNtHashArgMin) {
  auto minimizers = tb::NtHashArgMinMinimize(args_);
  auto precomputed_minimizers =
      tb::NtHashPrecomputedArgMinUnrolledMinimize(args_);

  EXPECT_EQ(minimizers, precomputed_minimizers);
}

TEST_F(MinimizeTest, NtHashArgMinRecoveryVsNtHashArgMin) {
  auto argmin_minimizers = tb::NtHashArgMinMinimize(args_);
  auto recovery_minimizers = tb::NtHashArgMinRecoveryMinimize(args_);

  EXPECT_EQ(argmin_minimizers, recovery_minimizers);
}

TEST_F(MinimizeTest, NtHashPrecomputedArgMinRecoveryVsNtHashArgMin) {
  auto argmin_minimizers = tb::NtHashArgMinMinimize(args_);
  auto recovery_minimizers = tb::NtHashPrecomputedArgMinRecoveryMinimize(args_);

  EXPECT_EQ(argmin_minimizers, recovery_minimizers);
}

TEST_F(MinimizeTest, NtHashPrecomputedArgMinUnrolledRecoveryVsNtHashArgMin) {
  auto argmin_minimizers = tb::NtHashArgMinMinimize(args_);
  auto recovery_minimizers =
      tb::NtHashPrecomputedArgMinUnrolledRecoveryMinimize(args_);

  EXPECT_EQ(argmin_minimizers, recovery_minimizers);
}

TEST_F(MinimizeTest, ArgMinRollingVsArgMin) {
  auto argmin_minimizers = tb::ArgMinMinimize(args_);
  auto rolling_minimizers = tb::ArgMinRollingMinimize(args_);

  EXPECT_EQ(argmin_minimizers, rolling_minimizers);
}

TEST_F(MinimizeTest, SplitWindowVsArgMin) {
  auto argmin_minimizers = tb::ArgMinMinimize(args_);
  auto split_minimizers = tb::SplitWindowMinimize(args_);

  EXPECT_EQ(argmin_minimizers, split_minimizers);
}
