cmake_minimum_required(VERSION 3.25)
project(cppbox LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 26)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

option(profile "Compile with profiling flags" OFF)

include(FetchContent)
FetchContent_Declare(
  doctest
  GIT_REPOSITORY https://github.com/doctest/doctest
  GIT_TAG v2.4.11)

FetchContent_Declare(
  nanobench
  GIT_REPOSITORY https://github.com/martinus/nanobench
  GIT_TAG v4.3.11)

FetchContent_MakeAvailable(doctest nanobench)

add_executable(cppbox-bench src/main.cc)
target_link_libraries(cppbox-bench PRIVATE doctest nanobench)
target_compile_options(
  cppbox-bench
  PRIVATE $<$<CONFIG:Debug,RelWithDebInfo>:-fsanitize=address
          -fsanitize=undefined -fno-omit-frame-pointer>
          $<$<CONFIG:Release>:-march=native>)
target_link_options(
  cppbox-bench PRIVATE $<$<CONFIG:Debug,RelWithDebInfo>:-fsanitize=address
  -fsanitize=undefined>)

if(profile)
  target_compile_options(
    cppbox-bench
    PRIVATE -g -fno-omit-frame-pointer
            $<$<CXX_COMPILER_ID:Clang>:-Rpass-analysis=loop-vectorize>
            $<$<CXX_COMPILER_ID:Clang>:-Rpass-analysis=loop-unroll>)
endif()
