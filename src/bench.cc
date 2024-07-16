#include <benchmark/benchmark.h>

#include "tb/algo.hpp"

namespace {

constexpr int kSeed = 42;

template <auto MinimizeFn> void BM_MinimizeW5K15(benchmark::State &state) {
  tb::MockSequence seq(state.range(0), kSeed);
  for (auto _ : state) {
    auto kmers = MinimizeFn({
        .seq = seq,
        .window_length = 5,
        .kmer_length = 15,
    });
    benchmark::DoNotOptimize(kmers.data());
  }
}

BENCHMARK_TEMPLATE(BM_MinimizeW5K15, tb::NaiveMinimize)
    ->RangeMultiplier(2)
    ->Range(1 << 10, 1 << 16);

BENCHMARK_TEMPLATE(BM_MinimizeW5K15, tb::DequeMinimize)
    ->RangeMultiplier(2)
    ->Range(1 << 10, 1 << 16);

BENCHMARK_TEMPLATE(BM_MinimizeW5K15, tb::InplaceMinimize)
    ->RangeMultiplier(2)
    ->Range(1 << 10, 1 << 16);

} // namespace
