//
// Created by Tiancheng Xie on 3/12/2021.
//

#ifndef DKG_VSS_FFT_H
#define DKG_VSS_FFT_H


#include <chrono>
#include "bn.h"
#include "gmpxx.h"
#include "gmp.h"

mpz_class* __dst[3];
bn::Ec1* __dst_ec[3];
mpz_class* twiddle_factor;
bn::Ec1* twiddle_factor_ec1;

void init_scratch_pad(int order)
{
    __dst[0] = new mpz_class[order];//(prime_field::field_element*)malloc(order * sizeof(prime_field::field_element));
    __dst[1] = new mpz_class[order];//(prime_field::field_element*)malloc(order * sizeof(prime_field::field_element));
    __dst[2] = new mpz_class[order];
    __dst_ec[0] = new bn::Ec1[order];
    __dst_ec[1] = new bn::Ec1[order];
    __dst_ec[2] = new bn::Ec1[order];
    twiddle_factor = new mpz_class[order];
    twiddle_factor_ec1 = new bn::Ec1[order];
}

void fast_fourier_transform(const mpz_class *coefficients, int coef_len, int order, mpz_class root_of_unity, mpz_class *result)
{
    mpz_class rot_mul[62];
    //note: malloc and free will not call the constructor and destructor, not recommended unless for efficiency
    //In sake of both memory and time efficiency, use the non-recursive version
    int lg_order = -1;
    rot_mul[0] = root_of_unity;
    for(int i = 0; i < 62; ++i)
    {
        if(i > 0)
            rot_mul[i] = rot_mul[i - 1] * rot_mul[i - 1] % p;
        if((1LL << i) == order)
        {
            lg_order = i;
        }
    }
    int lg_coef = -1;
    for(int i = 0; i < 62; ++i)
    {
        if((1LL << i) == coef_len)
        {
            lg_coef = i;
        }
    }
    assert(lg_order != -1 && lg_coef != -1);

    //we can merge both cases, but I just don't want to do so since it's easy to make mistake
    if(lg_coef > lg_order)
    {
        assert(false);
    }
    else
    {
        //initialize leaves
        int blk_sz = (order / coef_len);
        for(int j = 0; j < blk_sz; ++j)
        {
            for(int i = 0; i < coef_len; ++i)
            {
                __dst[lg_coef & 1][(j << lg_coef) | i] = coefficients[i];
            }
        }

        mpz_class *x_arr = new mpz_class[1 << lg_order];
        {
            for(int dep = lg_coef - 1; dep >= 0; --dep)
            {
                int blk_size = 1 << (lg_order - dep);
                int half_blk_size = blk_size >> 1;
                int cur = dep & 1;
                int pre = cur ^ 1;

                mpz_class x = 1;
                x_arr[0] = 1;
                for(int j = 1; j < blk_size; ++j)
                    x_arr[j] = x_arr[j - 1] * rot_mul[dep] % p;
                for(int k = 0; k < blk_size / 2; ++k)
                {
                    int double_k = (k) & (half_blk_size - 1);
                    for(int j = 0; j < (1 << dep); ++j)
                    {
                        mpz_class l_value = __dst[pre][double_k << (dep + 1) | j], r_value = x_arr[k] * __dst[pre][double_k << (dep + 1) | (1 << dep) | j];
                        __dst[cur][k << dep | j] = (l_value + r_value) % p;
                        __dst[cur][(k + blk_size / 2) << dep | j] = ((l_value - r_value) % p + p) % p;
                    }
                }
            }
        }
        delete[] x_arr;
    }

    for(int i = 0; i < order; ++i)
        result[i] = __dst[0][i];
}

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

void fast_fourier_transform(const bn::Ec1 *coefficients, int coef_len, int order, mpz_class root_of_unity, bn::Ec1 *result)
{
    mpz_class rot_mul[62];
    //note: malloc and free will not call the constructor and destructor, not recommended unless for efficiency
    //In sake of both memory and time efficiency, use the non-recursive version
    int lg_order = -1;
    rot_mul[0] = root_of_unity;
    for(int i = 0; i < 62; ++i)
    {
        if(i > 0)
            rot_mul[i] = rot_mul[i - 1] * rot_mul[i - 1] % p;
        if((1LL << i) == order)
        {
            lg_order = i;
        }
    }
    int lg_coef = -1;
    for(int i = 0; i < 62; ++i)
    {
        if((1LL << i) == coef_len)
        {
            lg_coef = i;
        }
    }
    assert(lg_order != -1 && lg_coef != -1);

    //we can merge both cases, but I just don't want to do so since it's easy to make mistake
    if(lg_coef > lg_order)
    {
        assert(false);
    }
    else
    {
        //initialize leaves
        int blk_sz = (order / coef_len);
        for(int j = 0; j < blk_sz; ++j)
        {
            for(int i = 0; i < coef_len; ++i)
            {
                __dst_ec[lg_coef & 1][(j << lg_coef) | i] = coefficients[i];
            }
        }

        mpz_class *x_arr = new mpz_class[1 << lg_order];
        {
            for(int dep = lg_coef - 1; dep >= 0; --dep)
            {
                int blk_size = 1 << (lg_order - dep);
                int half_blk_size = blk_size >> 1;
                int cur = dep & 1;
                int pre = cur ^ 1;

                mpz_class x = 1;
                x_arr[0] = 1;
                for(int j = 1; j < blk_size; ++j)
                    x_arr[j] = x_arr[j - 1] * rot_mul[dep] % p;
                for(int k = 0; k < blk_size / 2; ++k)
                {
                    int double_k = (k) & (half_blk_size - 1);
                    for(int j = 0; j < (1 << dep); ++j)
                    {
                        bn::Ec1 l_value = __dst_ec[pre][double_k << (dep + 1) | j], r_value = fastmul(__dst_ec[pre][double_k << (dep + 1) | (1 << dep) | j], x_arr[k]);
                        __dst_ec[cur][k << dep | j] = (l_value + r_value);
                        __dst_ec[cur][(k + blk_size / 2) << dep | j] = l_value - r_value;
                    }
                }
            }
        }
        delete[] x_arr;
    }

    for(int i = 0; i < order; ++i)
        result[i] = __dst_ec[0][i];
}

void inverse_fast_fourier_transform(bn::Ec1 *evaluations, int coef_len, int order, mpz_class root_of_unity, bn::Ec1 *dst)
{
    if(coef_len > order)
    {
        //more coefficient than evaluation
        fprintf(stderr, "Warning, Request do inverse fft with inefficient number of evaluations.");
        fprintf(stderr, "Will construct a polynomial with less order than required.");
        coef_len = order;
    }

    //assume coef_len <= order

    //subsample evalutions

    bn::Ec1 *sub_eval;
    bool need_free = false;
    if(coef_len != order)
    {
        need_free = true;
        sub_eval = new bn::Ec1[coef_len];
        for(int i = 0; i < coef_len; ++i)
        {
            sub_eval[i] = evaluations[i * (order / coef_len)];
        }
    }
    else
        sub_eval = evaluations;

    mpz_class new_rou = 1;
    for(int i = 0; i < order / coef_len; ++i)
        new_rou = new_rou * root_of_unity % p;
    order = coef_len;

    mpz_class inv_rou = 1, tmp = new_rou;
    int lg_order = -1;
    for(int i = 0; i < 62; ++i)
    {
        if((1LL << i) == order)
        {
            lg_order = i;
        }
    }
    int lg_coef = -1;
    for(int i = 0; i < 62; ++i)
    {
        if((1LL << i) == coef_len)
        {
            lg_coef = i;
        }
    }
    assert(lg_order != -1 && lg_coef != -1);

    for(int i = 0; i < lg_order; ++i)
    {
        inv_rou = inv_rou * tmp % p;
        tmp = tmp * tmp % p;
    }

    fast_fourier_transform(sub_eval, order, coef_len, inv_rou, dst);

    if(need_free)
        free(sub_eval);

    mpz_class inv_n = fastpow<mpz_class>(mpz_class(order), p - 2);

    for(int i = 0; i < coef_len; ++i)
    {
        dst[i] = dst[i] * inv_n;
    }
}


void fft_eval(const std::vector<mpz_class> &prover_coef, const mpz_class &omega, const int lg_players, std::vector<mpz_class> &eval_result)
{
    eval_result.resize(1 << lg_players);
    fast_fourier_transform(prover_coef.data(), prover_coef.size(), 1 << lg_players, omega, eval_result.data());
}

extern std::vector<bn::Ec1> B_rev_eval_ec;

void conv(const std::vector<mpz_class> &A, const std::vector<bn::Ec1> &B, std::vector<bn::Ec1> &ans)
{
    assert(A.size() == B.size());
    int d = A.size();
    ans.resize(2 * d);
    mpz_class rou = fastpow(root_of_unity, (1LL << (max_lg_order)) / (d * 2));
    mpz_class *eval_A;
    printf("conv start\n");
    //bn::Ec1 *eval_B, *eval_C;
    bn::Ec1 *eval_C;
    eval_A = new mpz_class[2 * d];
    //eval_B = new bn::Ec1[2 * d];
    eval_C = new bn::Ec1[2 * d];
    fast_fourier_transform(A.data(), d, 2 * d, rou, eval_A);
    //fast_fourier_transform(B.data(), d, 2 * d, rou, eval_B);
    for(int i = 0; i < 2 * d; ++i)
        eval_C[i] = B_rev_eval_ec[i] * eval_A[i];
    inverse_fast_fourier_transform(eval_C, 2 * d, 2 * d, rou, ans.data());
    printf("conv end\n");
}

#endif //DKG_VSS_FFT_H
