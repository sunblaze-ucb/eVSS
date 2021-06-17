//
// Created by Tiancheng Xie on 3/13/2021.
//

#include "bn.h"
#include "gmpxx.h"
#include "gmp.h"
#include "enviroment.h"
#include <chrono>

bn::Ec1 fastmul(bn::Ec1 x, mpz_class y)
{
    const int lg_windows_sz = 4;
    const int window_sz_int = 1 << lg_windows_sz;
    mpz_class window_sz = window_sz_int;
    bn::Ec1 x_arr[window_sz_int];
    x_arr[0] = g1 * 0;
    for(int i = 1; i < window_sz_int; ++i)
        x_arr[i] = x_arr[i - 1] + x;
    bn::Ec1 res = g1 * 0;
    static int y_bits[128];
    memset(y_bits, 0, sizeof(y_bits));
    int max_bit_pairs = 0;

    mpz_class y_parts_mpz[4] = {y, y >> 64, y >> 128, y >> 192};
    unsigned long long y_parts[4] = {y_parts_mpz[0].get_ui(), y_parts_mpz[1].get_ui(), y_parts_mpz[2].get_ui(), y_parts_mpz[3].get_ui()};


    for(int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 64 / lg_windows_sz; ++j) {
            y_bits[max_bit_pairs++] = y_parts[i] & (window_sz_int - 1);
            y_parts[i] /= window_sz_int;
        }
    }
    for(int i = max_bit_pairs - 1; i >= 0; --i)
    {
        bn::Ec1::dbl(res, res);
        bn::Ec1::dbl(res, res);
        bn::Ec1::dbl(res, res);
        bn::Ec1::dbl(res, res);
        res = res + x_arr[y_bits[i]];
    }
    return res;
}

const int maxn = 10000;

mpz_class rand_coef[maxn];
mie::Vuint rand_coef_vuint[maxn];

int main()
{
    init_enviroment();

    bn::Ec1 a = g1 * mie::Vuint(rand());
    bn::Ec1 res1 = g1 * 0;
    bn::Ec1 res2 = g1 * 0;

    for(int i = 0; i < maxn; ++i)
    {
        mpz_urandomm(rand_coef[i].get_mpz_t(), r_state, p.get_mpz_t());
    }

    auto t0 = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < maxn; ++i)
    {
        res1 = res1 + fastmul(a, rand_coef[i]);
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    printf("res1 generation time: %f\n", time_span.count());

    t0 = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < maxn; ++i)
    {
        res2 = res2 + a * rand_coef[i];
    }
    t1 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    printf("res2 generation time: %f\n", time_span.count());

    if(res1 == res2)
        printf("Pass\n");
    else
        printf("No Pass\n");

    return 0;
}
