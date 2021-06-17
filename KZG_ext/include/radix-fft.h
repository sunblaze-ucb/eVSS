//
// Created by Tiancheng Xie on 3/12/2021.
//

#ifndef DKG_VSS_FFT_H
#define DKG_VSS_FFT_H


#include <chrono>
#include "bn.h"
#include "gmpxx.h"
#include "gmp.h"
#include "enviroment.h"

mpz_class* __dst[3];
bn::Ec1* __dst_ec[3];
mpz_class* radix_fft_coef[22];
bn::Ec1* radix_fft_coef_ec[22];
mpz_class* radix_fft_result_U[22];
mpz_class* radix_fft_result_Z[22];
mpz_class* radix_fft_result_Z_[22];

bn::Ec1* radix_fft_result_U_ec[22];
bn::Ec1* radix_fft_result_Z_ec[22];
bn::Ec1* radix_fft_result_Z__ec[22];

mpz_class *twiddle_factor;

void init_scratch_pad(int order, mpz_class root_of_unity)
{
    __dst[0] = new mpz_class[order];//(prime_field::field_element*)malloc(order * sizeof(prime_field::field_element));
    __dst[1] = new mpz_class[order];//(prime_field::field_element*)malloc(order * sizeof(prime_field::field_element));
    __dst[2] = new mpz_class[order];
    __dst_ec[0] = new bn::Ec1[order];
    __dst_ec[1] = new bn::Ec1[order];
    __dst_ec[2] = new bn::Ec1[order];
    twiddle_factor = new mpz_class[order];

    mpz_class x = 1;
    for(int i = 0; i < order; ++i)
    {
        twiddle_factor[i] = x;
        x = x * root_of_unity % p;
    }
    for(int i = 0; (order >> i) > 0; ++i)
    {
        radix_fft_coef[i] = new mpz_class[order >> i];
        radix_fft_result_U[i] = new mpz_class[order >> i];
        radix_fft_result_Z[i] = new mpz_class[order >> i];
        radix_fft_result_Z_[i] = new mpz_class[order >> i];

        radix_fft_coef_ec[i] = new bn::Ec1[order >> i];
        radix_fft_result_U_ec[i] = new bn::Ec1[order >> i];
        radix_fft_result_Z_ec[i] = new bn::Ec1[order >> i];
        radix_fft_result_Z__ec[i] = new bn::Ec1[order >> i];

    }
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
    if(y == 1)
        return x;
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
int fft_cnt = 0;
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
                        fft_cnt++;
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

void conv(const std::vector<mpz_class> &A, const std::vector<bn::Ec1> &B, std::vector<bn::Ec1> &ans)
{
    assert(A.size() == B.size());
    int d = A.size();
    ans.resize(2 * d);
    mpz_class rou = fastpow(root_of_unity, (1LL << (max_lg_order)) / (d * 2));
    mpz_class *eval_A;
    printf("conv start\n");
    bn::Ec1 *eval_B, *eval_C;
    eval_A = new mpz_class[2 * d];
    eval_B = new bn::Ec1[2 * d];
    eval_C = new bn::Ec1[2 * d];
    fast_fourier_transform(A.data(), d, 2 * d, rou, eval_A);
    fast_fourier_transform(B.data(), d, 2 * d, rou, eval_B);
    for(int i = 0; i < 2 * d; ++i)
        eval_C[i] = eval_B[i] * eval_A[i];
    inverse_fast_fourier_transform(eval_C, 2 * d, 2 * d, rou, ans.data());
    printf("conv end\n");
}
#include <iostream>

inline mpz_class fastmod(mpz_class x)
{
    if(x < 0)
        return x + p;
    if(x >= p)
        return x - p;
    else
        return x;
}

void radix_fft(const mpz_class *coefficients, int coef_len, int order, mpz_class root_of_unity, mpz_class *result, int depth = 0)
{
    if(coef_len < 8)
    {
        fast_fourier_transform(coefficients, coef_len, order, root_of_unity, result);
    }
    else
    {
        for(int i = 0; i < order / 2; ++i)
        {
            radix_fft_coef[depth][i] = coefficients[i * 2];
        }
        radix_fft(radix_fft_coef[depth], coef_len / 2, order / 2, root_of_unity * root_of_unity % p, radix_fft_result_U[depth], depth + 1);
        for(int i = 0; i < order / 4; ++i)
        {
            radix_fft_coef[depth][i] = coefficients[i * 4 + 1];
        }
        radix_fft(radix_fft_coef[depth], coef_len / 4, order / 4, root_of_unity * root_of_unity * root_of_unity * root_of_unity % p, radix_fft_result_Z[depth], depth + 2);
        for(int i = 0; i < order / 4; ++i)
        {
            radix_fft_coef[depth][i] = coefficients[i * 4 + 3];
        }
        radix_fft(radix_fft_coef[depth], coef_len / 4, order / 4, root_of_unity * root_of_unity * root_of_unity * root_of_unity % p, radix_fft_result_Z_[depth], depth + 2);

        mpz_class rou_n_4 = fastpow<mpz_class>(root_of_unity, order / 4); // equals to -i
        for(int i = 0; i < order / 4; ++i)
        {
            const mpz_class &x = twiddle_factor[i << depth];
            const mpz_class &x3 = twiddle_factor[(3 * i) << depth];
            const mpz_class A = x * radix_fft_result_Z[depth][i] % p;
            const mpz_class B = x3 * radix_fft_result_Z_[depth][i] % p;
            const mpz_class ApB = fastmod(A + B), AmB = fastmod(A - B);
            mpz_class C = (rou_n_4 * AmB) % p;

            result[i] = fastmod(radix_fft_result_U[depth][i] + ApB);
            result[i + order / 2] = fastmod(radix_fft_result_U[depth][i] - ApB);
            result[i + order / 4] = fastmod(radix_fft_result_U[depth][i + order / 4] + C);
            result[i + 3 * order / 4] = fastmod(radix_fft_result_U[depth][i + order / 4] - C);
        }
    }
}
int cnt = 0;
void radix_fft(const bn::Ec1 *coefficients, int coef_len, int order, mpz_class root_of_unity, bn::Ec1 *result, int depth = 0)
{
    if(coef_len < 4)
    {
/*
        mpz_class rou = 1;
        for(int i = 0; i < order; ++i)
        {
            mpz_class x = rou;
            result[i] = coefficients[0];
            for(int j = 1; j < coef_len; ++j)
            {
                result[i] = result[i] + fastmul(coefficients[j], x);
                x = x * rou % p;
                cnt++;
            }
            rou = rou * root_of_unity;
            cnt++;
        }

        fast_fourier_transform(coefficients, coef_len, order, root_of_unity, result);
*/
        if(coef_len == 1)
        {
            for(int i = 0; i < order; ++i)
                result[i] = coefficients[0];
        }
        else
        {
            if(order == coef_len)
            {
                result[0] = coefficients[0] + coefficients[1];
                result[1] = coefficients[0] + coefficients[1] * root_of_unity;
            }
            else {
                mpz_class rou = 1;
                for (int i = 0; i < order; ++i) {
                    result[i] = coefficients[0] + coefficients[1] * rou;
                    rou = rou * root_of_unity;
                }
            }
        }
    }
    else
    {
        for(int i = 0; i < order / 2; ++i)
        {
            radix_fft_coef_ec[depth][i] = coefficients[i * 2];
        }
        radix_fft(radix_fft_coef_ec[depth], coef_len / 2, order / 2, root_of_unity * root_of_unity % p, radix_fft_result_U_ec[depth], depth + 1);
        for(int i = 0; i < order / 4; ++i)
        {
            radix_fft_coef_ec[depth][i] = coefficients[i * 4 + 1];
        }
        radix_fft(radix_fft_coef_ec[depth], coef_len / 4, order / 4, root_of_unity * root_of_unity * root_of_unity * root_of_unity % p, radix_fft_result_Z_ec[depth], depth + 2);
        for(int i = 0; i < order / 4; ++i)
        {
            radix_fft_coef_ec[depth][i] = coefficients[i * 4 + 3];
        }
        radix_fft(radix_fft_coef_ec[depth], coef_len / 4, order / 4, root_of_unity * root_of_unity * root_of_unity * root_of_unity % p, radix_fft_result_Z__ec[depth], depth + 2);

        mpz_class rou_n_4 = fastpow<mpz_class>(root_of_unity, order / 4);
        std::cout << rou_n_4 << std::endl;
        for(int i = 0; i < order / 4; ++i)
        {
            const mpz_class &x = twiddle_factor[i << depth];
            const mpz_class &x3 = twiddle_factor[(3 * i) << depth];
            const auto A = fastmul(radix_fft_result_Z_ec[depth][i], x);
            cnt++;
            const auto B = fastmul(radix_fft_result_Z__ec[depth][i], x3);
            cnt++;
            const auto ApB = A + B, AmB = A - B;

            const auto C = fastmul(AmB, rou_n_4);
            cnt++;

            result[i] = (radix_fft_result_U_ec[depth][i] + ApB);
            result[i + order / 2] = (radix_fft_result_U_ec[depth][i] - ApB);
            result[i + order / 4] = (radix_fft_result_U_ec[depth][i + order / 4] + C);
            result[i + 3 * order / 4] = (radix_fft_result_U_ec[depth][i + order / 4] - C);
        }
    }
}

#endif //DKG_VSS_FFT_H
