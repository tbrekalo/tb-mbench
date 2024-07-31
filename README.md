# Minimizer benchmarks
Different implementations for finding robust random minimizers in genomic data.

## Running benchmarks

### Dependencites
- C++20 compliant compiler
  - Support for [C++17 execution policy](https://en.cppreference.com/w/cpp/algorithm/execution_policy_tag_t)
- CMake
- git

- [googletest](https://github.com/google/googletest)
- [benchmark](https://github.com/google/benchmark)
- [eve](https://github.com/jfalcou/eve)

### Build
```bash
git clone https://github.com/tbrekalo/tb-mbench
cd tb-mbench

cmake -Bbuild -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

### Run
```bash
./build/bin/bench
```

## Results

### Benchmarking machine
- 13th Gen Intel(R) Core(TM) i5-13500
- locked at 2.5Ghz
- window length = 11
- kmer length = 21

### Figures
<matplotlib.legend.Legend at 0x7fdfecae4dd0>
| function                                            |   avg ns per bp |
|:----------------------------------------------------|----------------:|
| tb::NaiveMinimize                                   |       241.621   |
| tb::ArgMinDuffMinimize                              |        28.8162  |
| tb::RingMinimize                                    |        26.2489  |
| tb::ArgMinEveMinimize                               |        25.3979  |
| tb::NtHashArgMinMinimize                            |        23.6542  |
| tb::DequeMinimize                                   |        23.3451  |
| tb::NtHashPrecomputedArgMinMinimize                 |        23.2418  |
| tb::ArgMinMinimize                                  |        22.69    |
| tb::InplaceMinimize                                 |        20.8377  |
| tb::NtHashPrecomputedArgMinUnrolledMinimize         |        12.0209  |
| tb::SplitWindowMinimize                             |        11.8228  |
| tb::ArgMinUnrolledMinimize                          |        11.5503  |
| tb::NtHashArgMinRecoveryMinimize                    |        11.4736  |
| tb::ArgMinRecoveryEveMinimize                       |        10.9206  |
| tb::ArgMinRollingMinimize                           |        10.8241  |
| tb::ArgMinRecoveryMinimize                          |        10.587   |
| tb::NtHashPrecomputedArgMinRecoveryMinimize         |        10.5757  |
| tb::ArgMinRecoveryDuffMinimize                      |         9.14772 |
| tb::ArgMinRecoveryUnrolledMinimize                  |         8.77331 |
| tb::NtHashPrecomputedArgMinUnrolledRecoveryMinimize |         8.54623 |
![](misc/perf.png)

## Reference
- [Winnowing: Local Algorithms for Document Fingerprinting](http://dx.doi.org/10.1145/872769.872770)
- [Reducing storage requirements for biological sequence comparison](https://doi.org/10.1093/bioinformatics/bth408)
- [CURIOUS CODING NtHash](https://curiouscoding.nl/posts/nthash/)
- [Mohamadi, Hamid, Justin Chu, Benjamin P. Vandervalk, and Inanc Birol. 2016. “Nthash: Recursive Nucleotide Hashing.” Bioinformatics 32 (22): 3492–94.](http://dx.doi.org/10.1093/bioinformatics/btw397)
- [Codeforces split window blog post](https://codeforces.com/blog/entry/71687)