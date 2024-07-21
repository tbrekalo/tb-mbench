cmake_minimum_required(VERSION 3.25)
project(mbench LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 26)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

option(profile "Compile with profiling flags" OFF)

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest
  GIT_TAG v1.14.0)

FetchContent_Declare(
  benchmark
  GIT_REPOSITORY https://github.com/google/benchmark
  GIT_TAG v1.8.4)

FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
  FetchContent_Populate(googletest)
  add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR}
                   EXCLUDE_FROM_ALL)
endif()

FetchContent_GetProperties(benchmark)
if(NOT benchmark_POPULATED)
  FetchContent_Populate(benchmark)
  add_subdirectory(${benchmark_SOURCE_DIR} ${benchmark_BINARY_DIR}
                   EXCLUDE_FROM_ALL)
endif()

set(EVE_BUILD_TEST OFF)
set(EVE_BUILD_BENCHMARKS OFF)
set(EVE_BUILD_DOCUMENTATION OFF)

FetchContent_Declare(
  eve
  GIT_REPOSITORY https://github.com/jfalcou/eve
  GIT_TAG v2023.02.15)

FetchContent_GetProperties(eve)
if(NOT eve_POPULATED)
  FetchContent_Populate(eve)
  add_subdirectory(${eve_SOURCE_DIR} ${eve_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

add_library(lib src/algo.cc src/data.cc src/nthash.cc src/ring.cc)
target_link_libraries(lib eve::eve)
target_include_directories(lib PUBLIC include)
target_compile_options(
  lib
  PUBLIC $<$<CONFIG:Debug,RelWithDebInfo>:-fsanitize=address
         -fsanitize=undefined -fno-omit-frame-pointer>
         $<$<CONFIG:Release>:-march=native>)
target_link_options(
  lib PUBLIC $<$<CONFIG:Debug,RelWithDebInfo>:-fsanitize=address
  -fsanitize=undefined>)

if(profile)
  target_compile_options(
    lib
    PUBLIC -g -fno-omit-frame-pointer
           $<$<CXX_COMPILER_ID:Clang>:-Rpass-analysis=loop-vectorize>
           $<$<CXX_COMPILER_ID:Clang>:-Rpass-analysis=loop-unroll>)
endif()

add_executable(bench src/bench.cc)
target_link_libraries(bench PRIVATE benchmark::benchmark_main lib)

add_executable(test src/test.cc)
target_link_libraries(test PRIVATE GTest::gtest_main lib)
