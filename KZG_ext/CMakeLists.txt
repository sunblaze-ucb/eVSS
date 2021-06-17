project(DKG_VSS)
cmake_minimum_required(VERSION 3.16)

add_compile_definitions(BN_SUPPORT_SNARK MIE_ATE_USE_GMP)
set(CMAKE_CXX_FLAGS "-march=native -std=c++11 -lstdc++ -O3")
include_directories(include include/ate-pairing/include)
link_directories(include/ate-pairing/lib)

add_executable(KZG_std src/kzg_main.cpp include/algebra/polynomial.h include/KZG_verifier.h include/test_point.hpp include/trusted_setup.h include/enviroment.h include/KZG_std_prover.h include/KZG_improved_prover.h include/fft.h include/KZG_commit.h)

add_executable(benchmark src/benchmark_algebra.cpp include/algebra/polynomial.h include/KZG_verifier.h include/test_point.hpp include/trusted_setup.h include/enviroment.h include/KZG_std_prover.h include/KZG_improved_prover.h include/fft.h include/KZG_commit.h)

add_executable(radix_fft_test src/radix_fft_test.cpp)
target_link_libraries(KZG_std zm.a gmpxx gmp)
target_link_libraries(benchmark zm.a gmpxx gmp)
target_link_libraries(radix_fft_test zm.a gmpxx gmp)