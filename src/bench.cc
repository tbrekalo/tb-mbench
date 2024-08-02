#include <benchmark/benchmark.h>

#include "tb/algo.hpp"

namespace {

constexpr int kSeed = 42;
constexpr std::size_t kNBasesSmall = 1'000uz;
constexpr std::size_t kNBasesLarge = 1'000'000uz;

template <auto MinimizeFn>
void BM_Minimize(benchmark::State& state) {
  for (auto _ : state) {
    state.PauseTiming();
    tb::MockSequence seq(state.range(0), kSeed);
    state.ResumeTiming();

    auto kmers = MinimizeFn({
        .seq = seq,
        .window_length = 11,
        .kmer_length = 21,
    });

    benchmark::DoNotOptimize(kmers.data());
  }
}

std::vector<std::vector<std::int64_t>> kArgList = {{kNBasesLarge}};

// Reference
BENCHMARK_TEMPLATE(BM_Minimize, tb::NaiveMinimize)->ArgsProduct(kArgList);

// Deque based
BENCHMARK_TEMPLATE(BM_Minimize, tb::DequeMinimize)->ArgsProduct(kArgList);

// Arg min based
BENCHMARK_TEMPLATE(BM_Minimize, tb::ArgMinMinimize)->ArgsProduct(kArgList);
BENCHMARK_TEMPLATE(BM_Minimize, tb::ArgMinUnrolledMinimize)
    ->ArgsProduct(kArgList);
BENCHMARK_TEMPLATE(BM_Minimize, tb::NtHashPrecomputedArgMinMinimize)
    ->ArgsProduct(kArgList);

// Arg min recovery based
BENCHMARK_TEMPLATE(BM_Minimize, tb::ArgMinRecoveryMinimize)
    ->ArgsProduct(kArgList);
BENCHMARK_TEMPLATE(BM_Minimize, tb::ArgMinRecoveryUnrolledMinimize)
    ->ArgsProduct(kArgList);
BENCHMARK_TEMPLATE(BM_Minimize,
                   tb::NtHashPrecomputedArgMinUnrolledRecoveryMinimize)
    ->ArgsProduct(kArgList);

// Split window
BENCHMARK_TEMPLATE(BM_Minimize, tb::SplitWindowMinimize)->ArgsProduct(kArgList);

}  // namespace
