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
        x = x * omega % p;
    }
    printf("proof generate done\n");
    return ret;
}

#endif //DKG_VSS_KZG_STD_PROVER_H
