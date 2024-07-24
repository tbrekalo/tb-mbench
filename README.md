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
| tb::NaiveMinimize                           |       241.547   |
| tb::NtHashArgMinMinimize                    |        26.454   |
| tb::RingMinimize                            |        26.2809  |
| tb::DequeMinimize                           |        23.495   |
| tb::NtHashPrecomputedArgMinMinimize         |        22.7503  |
| tb::InplaceMinimize                         |        21.313   |
| tb::NtHashArgMinRecoveryMinimize            |        20.6429  |
| tb::NtHashPrecomputedArgMinRecoveryMinimize |        16.5655  |
| tb::ArgMinMinimize                          |        15.9429  |
| tb::ArgMinRecoveryEveMinimize               |        12.8927  |
| tb::ArgMinRecoveryMinimize                  |         9.81566 |
![](misc/perf.png)

## Reference
- [Winnowing: Local Algorithms for Document Fingerprinting](http://dx.doi.org/10.1145/872769.872770)
- [Reducing storage requirements for biological sequence comparison](https://doi.org/10.1093/bioinformatics/bth408)
- [CURIOUS CODING NtHash](https://curiouscoding.nl/posts/nthash/)
- [Mohamadi, Hamid, Justin Chu, Benjamin P. Vandervalk, and Inanc Birol. 2016. “Nthash: Recursive Nucleotide Hashing.” Bioinformatics 32 (22): 3492–94.](http://dx.doi.org/10.1093/bioinformatics/btw397)