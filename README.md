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
| function                                    |   avg ns per bp |
|:--------------------------------------------|----------------:|
| tb::NaiveMinimize                           |       241.791   |
| tb::RingMinimize                            |        26.3036  |
| tb::DequeMinimize                           |        23.4914  |
| tb::InplaceMinimize                         |        21.3016  |
| tb::NtHashArgMinMinimize                    |        15.6999  |
| tb::ArgMinMinimize                          |        15.0884  |
| tb::NtHashPrecomputedArgMinMinimize         |        15.0337  |
| tb::ArgMinRecoveryEveMinimize               |        11.913   |
| tb::NtHashArgMinRecoveryMinimize            |         9.45148 |
| tb::ArgMinRecoveryMinimize                  |         8.87511 |
| tb::NtHashPrecomputedArgMinRecoveryMinimize |         8.81211 |
![](misc/perf.png)

## Reference
- [Winnowing: Local Algorithms for Document Fingerprinting](http://dx.doi.org/10.1145/872769.872770)
- [Reducing storage requirements for biological sequence comparison](https://doi.org/10.1093/bioinformatics/bth408)
- [CURIOUS CODING NtHash](https://curiouscoding.nl/posts/nthash/)
- [Mohamadi, Hamid, Justin Chu, Benjamin P. Vandervalk, and Inanc Birol. 2016. “Nthash: Recursive Nucleotide Hashing.” Bioinformatics 32 (22): 3492–94.](http://dx.doi.org/10.1093/bioinformatics/btw397)