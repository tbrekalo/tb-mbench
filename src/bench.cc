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

BENCHMARK_TEMPLATE(BM_Minimize, tb::NaiveMinimize)->ArgsProduct(kArgList);
BENCHMARK_TEMPLATE(BM_Minimize, tb::DequeMinimize)->ArgsProduct(kArgList);
BENCHMARK_TEMPLATE(BM_Minimize, tb::InplaceMinimize)->ArgsProduct(kArgList);
BENCHMARK_TEMPLATE(BM_Minimize, tb::RingMinimize)->ArgsProduct(kArgList);
BENCHMARK_TEMPLATE(BM_Minimize, tb::ArgMinMinimize)->ArgsProduct(kArgList);
BENCHMARK_TEMPLATE(BM_Minimize, tb::NtHashArgMinMinimize)
    ->ArgsProduct(kArgList);
BENCHMARK_TEMPLATE(BM_Minimize, tb::NtHashPrecomputedArgMinMinimize)
    ->ArgsProduct(kArgList);
BENCHMARK_TEMPLATE(BM_Minimize, tb::ArgMinRecoveryMinimize)
    ->ArgsProduct(kArgList);
BENCHMARK_TEMPLATE(BM_Minimize, tb::ArgMinRecoveryEveMinimize)
    ->ArgsProduct(kArgList);
BENCHMARK_TEMPLATE(BM_Minimize, tb::NtHashArgMinRecoveryMinimize)
    ->ArgsProduct(kArgList);
BENCHMARK_TEMPLATE(BM_Minimize, tb::NtHashPrecomputedArgMinRecoveryMinimize)
    ->ArgsProduct(kArgList);

}  // namespace
