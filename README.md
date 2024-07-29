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
| function                                            |   avg ns per bp |
|:----------------------------------------------------|----------------:|
| tb::NaiveMinimize                                   |       241.694   |
| tb::ArgMinDuffMinimize                              |        28.4622  |
| tb::RingMinimize                                    |        26.2314  |
| tb::ArgMinEveMinimize                               |        25.3984  |
| tb::DequeMinimize                                   |        23.36    |
| tb::NtHashArgMinMinimize                            |        23.1313  |
| tb::ArgMinMinimize                                  |        22.5943  |
| tb::NtHashPrecomputedArgMinMinimize                 |        22.4447  |
| tb::InplaceMinimize                                 |        20.8296  |
| tb::SplitWindowMinimize                             |        18.5789  |
| tb::NtHashPrecomputedArgMinUnrolledMinimize         |        11.6239  |
| tb::ArgMinUnrolledMinimize                          |        11.4983  |
| tb::ArgMinRollingMinimize                           |        11.107   |
| tb::NtHashArgMinRecoveryMinimize                    |        11.0996  |
| tb::ArgMinRecoveryEveMinimize                       |        10.8627  |
| tb::NtHashPrecomputedArgMinRecoveryMinimize         |        10.4391  |
| tb::ArgMinRecoveryMinimize                          |        10.3098  |
| tb::NtHashPrecomputedArgMinUnrolledRecoveryMinimize |         8.9795  |
| tb::ArgMinRecoveryDuffMinimize                      |         8.7492  |
| tb::ArgMinRecoveryUnrolledMinimize                  |         8.74724 |
![](misc/perf.png)

## Reference
- [Winnowing: Local Algorithms for Document Fingerprinting](http://dx.doi.org/10.1145/872769.872770)
- [Reducing storage requirements for biological sequence comparison](https://doi.org/10.1093/bioinformatics/bth408)
- [CURIOUS CODING NtHash](https://curiouscoding.nl/posts/nthash/)
- [Mohamadi, Hamid, Justin Chu, Benjamin P. Vandervalk, and Inanc Birol. 2016. “Nthash: Recursive Nucleotide Hashing.” Bioinformatics 32 (22): 3492–94.](http://dx.doi.org/10.1093/bioinformatics/btw397)
- [Codeforces split window blog post](https://codeforces.com/blog/entry/71687)