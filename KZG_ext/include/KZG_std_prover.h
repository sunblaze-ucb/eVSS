//
// Created by Tiancheng Xie on 3/10/2021.
//

#ifndef DKG_VSS_KZG_STD_PROVER_H
#define DKG_VSS_KZG_STD_PROVER_H

#include <vector>
#include "bn.h"
#include <iostream>
#include "enviroment.h"

extern std::vector<mpz_class> prover_coef;

std::vector<single_proof> generate_proof(int lg_degree, int lg_players, const std::vector<bn::Ec1> &pp)
{
    std::vector<single_proof> ret;
    mpz_class x = 1;
    auto omega = fastpow(root_of_unity, 1LL << (max_lg_order - lg_players));
    //std::cout << "omega " << omega << " omega^2 " << omega * omega % p << std::endl;
    mpz_class *new_coef = new mpz_class [1 << lg_degree];
    new_coef[(1 << lg_degree) - 1] = 0;
    for(int i = 0; i < (1 << lg_players); ++i)
    {
        mpz_class eval_result = 0;
        mpz_class current_x = 1;
        for(int j = 0; j < (1 << lg_degree); ++j)
        {
            eval_result = (eval_result + prover_coef[j] * current_x) % p;
            current_x = current_x * x % p;
            //std::cout << "prover_coef[" << j << "] = " << prover_coef[j] << std::endl;
        }

        single_proof current_proof;
        current_proof.ans = eval_result;

        for(int j = (1 << lg_degree) - 2; j >= 0; --j)
        {
            new_coef[j] = new_coef[j + 1];
            new_coef[j] = new_coef[j] * x;
            new_coef[j] = (new_coef[j] + prover_coef[j + 1]) % p;
            //std::cout << "new_coef[" << j << "] = " << new_coef[j] << std::endl;
        }

        current_proof.pi = g1 * 0;
        for(int j = 0; j < (1 << lg_degree) - 1; ++j)
        {
            mie::Vuint coef_vuint(new_coef[j].get_str().c_str());
            current_proof.pi = current_proof.pi + pp[j] * coef_vuint;
        }
        ret.push_back(current_proof);
/*
        mpz_class rand_test;
        rand_test = tau;
        //std::cout << "rand_test " << rand_test << std::endl;
        //std::cout << "omega " << x << std::endl;
        mpz_class var = 1;
        mpz_class res1 = 0, res_tau = 0, res_s_i = 0;
        for(int j = 0; j < (1 << lg_degree) - 1; ++j)
        {
            res1 = (res1 + new_coef[j] * var) % p;
            var = (var * rand_test) % p;
        }
        var = 1;
        for(int j = 0; j < (1 << lg_degree); ++j)
        {
            res_tau = (res_tau + prover_coef[j] * var) % p;
            var = var * rand_test % p;
        }
        res_s_i = eval_result;

        mpz_class de = normalize(rand_test - x);
        //std::cout << de << std::endl;
        assert(de >= 0 && de < p);
        mpz_class res_div = normalize(res_tau - res_s_i) * fastpow<mpz_class>(de, p - 2) % p;
        //std::cout << "res1 via new_coef " << res1 << std::endl << "res_div via div " << res_div << std::endl;
        //std::cout << "res_tau " << res_tau << std::endl << "res_s_i " << res_s_i << " " << std::endl << "normalized " << normalize(res_tau - res_s_i) << std::endl;
        assert(res1 == res_div);

        assert(g1 * mie::Vuint(res1.get_str().c_str()) == current_proof.pi);


        {
            bn::Fp12 e1, e2;
            mie::Vuint vuint_y(eval_result.get_str().c_str());
            mie::Vuint vuint_x(x.get_str().c_str());
            auto e1x = g2;
            auto c = g1 * 0;
            mpz_class eval_tau = 0;
            mpz_class current_tau = 1;
            for(int j = 0; j < (1 << lg_degree); ++j)
            {
                mie::Vuint coef_vuint(prover_coef[j].get_str().c_str());
                c = c + pp[j] * coef_vuint;
                eval_tau = eval_tau + prover_coef[j] * current_tau;
                current_tau *= tau;
            }
            assert(g1 * mie::Vuint(eval_tau.get_str().c_str()) == c);
            bn::Ec1 e1y = c - (g1 * vuint_y);
            bn::Ec2 e2x = (g2 * mie::Vuint(tau.get_str().c_str())) - (g2 * vuint_x);
            auto e2y = current_proof.pi;

            bn::opt_atePairing(e1, e1x, e1y);
            bn::opt_atePairing(e2, e2x, e2y);

            assert(e1 == e2);
        }

*/
        x = x * omega % p;
    }
    printf("proof generate done\n");
    return ret;
}

#endif //DKG_VSS_KZG_STD_PROVER_H
