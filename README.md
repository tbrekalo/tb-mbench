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

### Figures
| function                |   avg ns per bp |
|:------------------------|----------------:|
| tb::NaiveMinimize       |         90.51   |
| tb::ArgminEveMinimize   |         31.9355 |
| tb::RingMinimize        |         22.4501 |
| tb::InplaceMinimize     |         19.7748 |
| tb::DequeMinimize       |         19.0838 |
| tb::ArgminUnseqMinimize |         17.7876 |
| tb::ArgminMinimize      |         17.773  |
![](misc/perf.png)

## Reference
- [Winnowing: Local Algorithms for Document Fingerprinting](http://dx.doi.org/10.1145/872769.872770)
- [Reducing storage requirements for biological sequence comparison](https://doi.org/10.1093/bioinformatics/bth408)