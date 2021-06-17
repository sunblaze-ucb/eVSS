//
// Created by Tiancheng Xie on 3/10/2021.
//

#ifndef DKG_VSS_TRUSTED_SETUP_H
#define DKG_VSS_TRUSTED_SETUP_H

#include <vector>
#include "bn.h"
#include "test_point.hpp"
#include "enviroment.h"
#include <utility>
#include <gmp.h>
#include <gmpxx.h>
#include "fft.h"
extern bn::Ec1 g1;
extern bn::Ec2 g2;
extern mpz_class p;
extern gmp_randstate_t r_state;

bn::Ec1 g_baby_step[256 / 15 + 1][1 << 15];

const int multi_scalar_w = 8;
class multi_scalar_state
{
public:
    bn::Ec1 value[(1 << multi_scalar_w) - 3];
};

void KeyGen_preprocessing(bn::Ec1 g)
{
    bn::Ec1 g_pow = g;
    for(int i = 0; i < 256 / 15 + 1; ++i)
    {
        g_baby_step[i][0] = g * 0;
        for(int j = 1; j < (1 << 15); ++j)
            g_baby_step[i][j] = g_baby_step[i][j - 1] + g_pow;
        g_pow = g_pow * (1 << 15);
    }
}

bn::Ec1 g1_exp(mpz_class a)
{
    bn::Ec1 ret = g1 * 0;
    int length = mpz_sizeinbase(a.get_mpz_t(), 2);
    for(int i = 0; i < length / 15 + 1; ++i)
    {
        mpz_class mask_mpz = (a >> (i * 15)) & ((1 << 15) - 1);
        int mask = mpz_get_si(mask_mpz.get_mpz_t());
        ret = ret + g_baby_step[i][mask];
    }
    return ret;
}

std::vector<bn::Ec1> B_rev_eval_ec;

std::pair<std::vector<bn::Ec1>, bn::Ec2 > trusted_setup_ec1(int d)
{
    mpz_class tau;
    KeyGen_preprocessing(g1);
    std::vector<bn::Ec1> tau_0;
    bn::Ec1 tau_0_x = g1;

    mpz_class x = 1;
    mpz_urandomm(tau.get_mpz_t(), r_state, p.get_mpz_t());
    tau = tau % p;
    //std::cout << "setup tau " << tau << std::endl;


    for(int i = 0; i < d; ++i)
    {
        tau_0_x = g1_exp(x);
        tau_0.push_back(tau_0_x);
        assert(g1_exp(x) == g1 * x);

        x = x * tau % p;
    }
    printf("Trusted setup done\n");

    std::vector<mpz_class> pp_rev(d), pp_rev_eval(2 * d);
    B_rev_eval_ec.resize(2 * d);

    x = 1;
    for(int i = d - 2; i >= 0; --i)
    {
        pp_rev[i] = x;
        x = x * tau % p;
    }
    pp_rev[d - 1] = 0;
    mpz_class rou = fastpow(root_of_unity, (1LL << (max_lg_order)) / (d * 2));
    fast_fourier_transform(pp_rev.data(), d, 2 * d, rou, pp_rev_eval.data());

    for(int i = 0; i < 2 * d; ++i)
    {
        B_rev_eval_ec[i] = g1 * pp_rev_eval[i];
    }


    return std::make_pair(tau_0, g2 * tau);
}

#endif //DKG_VSS_TRUSTED_SETUP_H
