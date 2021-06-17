//
// Created by Tiancheng Xie on 4/14/2021.
//

#include "radix_fft_test.h"
#include <cstdio>
#include "trusted_setup.h"
#include "radix-fft.h"
#include "enviroment.h"
#include <gmp.h>
#include <gmpxx.h>
#include <chrono>

bn::Ec1 *coef, *result, *result2;
extern int cnt;
extern int fft_cnt;
int main()
{
    int lg_N = 10;
    int N = 1 << lg_N;
    init_enviroment();
    auto omega = fastpow(root_of_unity, 1LL << (max_lg_order - lg_N));
    init_scratch_pad(N, omega);
    coef = new bn::Ec1[N];
    result = new bn::Ec1[N];
    result2 = new bn::Ec1[N];
    for(int i = 0; i < N; ++i)
    {
        coef[i] = g1 * rand();
    }
    auto t0 = std::chrono::high_resolution_clock::now();
    fast_fourier_transform(coef, N, N, omega, result);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    std::cout << "std fft time: " << time_span.count() << std::endl;

    t0 = std::chrono::high_resolution_clock::now();
    radix_fft(coef, N, N, omega, result2);
    t1 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    std::cout << "radix fft time: " << time_span.count() << std::endl;
    std::cout << cnt << " " << fft_cnt << std::endl;
    for(int i = 0; i < N; ++i)
    {
        assert(result[i] == result2[i]);
    }
    return 0;
}