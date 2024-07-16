# Minimizer benchmarks
Different implementations for finding robust random minimizers in genomic data.

## Running benchmarks

### Dependencites
- C++20 compliant compiler
- CMake
- git

- [googletest](https://github.com/google/googletest)
- [benchmark](https://github.com/google/benchmark)

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
| function            |   avg ns per bp |   std ns per bp |
|:--------------------|----------------:|----------------:|
| tb::NaiveMinimize   |         89.2503 |        0.765395 |
| tb::DequeMinimize   |         13.9247 |        4.51882  |
| tb::InplaceMinimize |         12.465  |        5.3303   |

![](misc/perf.png)

## Reference
- [Winnowing: Local Algorithms for Document Fingerprinting](http://dx.doi.org/10.1145/872769.872770)
- [Reducing storage requirements for biological sequence comparison](https://doi.org/10.1093/bioinformatics/bth408)