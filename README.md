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
| tb::NaiveMinimize                                   |       241.478   |
| tb::RingMinimize                                    |        26.2147  |
| tb::ArgMinEveMinimize                               |        25.3513  |
| tb::DequeMinimize                                   |        23.4849  |
| tb::InplaceMinimize                                 |        21.4213  |
| tb::NtHashArgMinMinimize                            |        15.5733  |
| tb::NtHashPrecomputedArgMinMinimize                 |        14.9972  |
| tb::ArgMinMinimize                                  |        14.952   |
| tb::NtHashPrecomputedArgMinUnrolledRecoveryMinimize |        13.5336  |
| tb::NtHashPrecomputedArgMinUnrolledMinimize         |        13.4937  |
| tb::ArgMinRecoveryUnrolledMinimize                  |        13.4047  |
| tb::ArgMinUnrolledMinimize                          |        13.357   |
| tb::ArgMinRecoveryEveMinimize                       |        11.8969  |
| tb::NtHashArgMinRecoveryMinimize                    |         9.56785 |
| tb::NtHashPrecomputedArgMinRecoveryMinimize         |         8.97545 |
| tb::ArgMinRecoveryMinimize                          |         8.87539 |
![](misc/perf.png)

## Reference
- [Winnowing: Local Algorithms for Document Fingerprinting](http://dx.doi.org/10.1145/872769.872770)
- [Reducing storage requirements for biological sequence comparison](https://doi.org/10.1093/bioinformatics/bth408)
- [CURIOUS CODING NtHash](https://curiouscoding.nl/posts/nthash/)
- [Mohamadi, Hamid, Justin Chu, Benjamin P. Vandervalk, and Inanc Birol. 2016. “Nthash: Recursive Nucleotide Hashing.” Bioinformatics 32 (22): 3492–94.](http://dx.doi.org/10.1093/bioinformatics/btw397)