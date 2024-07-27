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
| tb::NaiveMinimize                                   |       241.667   |
| tb::RingMinimize                                    |        26.1466  |
| tb::ArgMinEveMinimize                               |        25.3506  |
| tb::DequeMinimize                                   |        23.5147  |
| tb::NtHashArgMinMinimize                            |        22.9285  |
| tb::NtHashPrecomputedArgMinMinimize                 |        22.3194  |
| tb::ArgMinMinimize                                  |        22.1474  |
| tb::InplaceMinimize                                 |        21.4186  |
| tb::NtHashPrecomputedArgMinUnrolledMinimize         |        11.6007  |
| tb::ArgMinUnrolledMinimize                          |        11.5099  |
| tb::NtHashArgMinRecoveryMinimize                    |        10.9382  |
| tb::ArgMinRecoveryEveMinimize                       |        10.5707  |
| tb::NtHashPrecomputedArgMinRecoveryMinimize         |        10.2953  |
| tb::ArgMinRecoveryMinimize                          |        10.2067  |
| tb::NtHashPrecomputedArgMinUnrolledRecoveryMinimize |         8.36242 |
| tb::ArgMinRecoveryUnrolledMinimize                  |         8.28662 |
![](misc/perf.png)

## Reference
- [Winnowing: Local Algorithms for Document Fingerprinting](http://dx.doi.org/10.1145/872769.872770)
- [Reducing storage requirements for biological sequence comparison](https://doi.org/10.1093/bioinformatics/bth408)
- [CURIOUS CODING NtHash](https://curiouscoding.nl/posts/nthash/)
- [Mohamadi, Hamid, Justin Chu, Benjamin P. Vandervalk, and Inanc Birol. 2016. “Nthash: Recursive Nucleotide Hashing.” Bioinformatics 32 (22): 3492–94.](http://dx.doi.org/10.1093/bioinformatics/btw397)