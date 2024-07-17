#include <benchmark/benchmark.h>

#include "tb/algo.hpp"

namespace {

constexpr int kSeed = 42;
constexpr std::size_t kNBases = 1'000'000uz;

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

BENCHMARK_TEMPLATE(BM_MinimizeW5K15, tb::NaiveMinimize)->Arg(kNBases);
BENCHMARK_TEMPLATE(BM_MinimizeW5K15, tb::DequeMinimize)->Arg(kNBases);
BENCHMARK_TEMPLATE(BM_MinimizeW5K15, tb::InplaceMinimize)->Arg(kNBases);
BENCHMARK_TEMPLATE(BM_MinimizeW5K15, tb::RingMinimize)->Arg(kNBases);
BENCHMARK_TEMPLATE(BM_MinimizeW5K15, tb::ArgminMinimize)->Arg(kNBases);
BENCHMARK_TEMPLATE(BM_MinimizeW5K15, tb::ArgminUnseqMinimize)->Arg(kNBases);
BENCHMARK_TEMPLATE(BM_MinimizeW5K15, tb::ArgminEveMinimize)->Arg(kNBases);

} // namespace
