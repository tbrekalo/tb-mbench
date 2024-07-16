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

BENCHMARK_TEMPLATE(BM_MinimizeW5K15, tb::NaiveMinimize)->Arg(1'000'000);
BENCHMARK_TEMPLATE(BM_MinimizeW5K15, tb::DequeMinimize)->Arg(1'000'000);
BENCHMARK_TEMPLATE(BM_MinimizeW5K15, tb::InplaceMinimize)->Arg(1'000'000);
BENCHMARK_TEMPLATE(BM_MinimizeW5K15, tb::RingMinimize)->Arg(1'000'000);

} // namespace
