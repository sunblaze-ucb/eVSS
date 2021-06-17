//
// Created by Tiancheng Xie on 3/12/2021.
//

#ifndef DKG_VSS_KZG_IMPROVED_PROVER_H
#define DKG_VSS_KZG_IMPROVED_PROVER_H
#include <vector>
#include "bn.h"
#include "gmpxx.h"
#include "gmp.h"
#include "enviroment.h"
#include "fft.h"

extern std::vector<mpz_class> prover_coef;

std::vector<single_proof> generate_proof(int lg_degree, int lg_players, const std::vector<bn::Ec1> &pp)
{
    std::vector<single_proof> ret;
    auto omega = fastpow(root_of_unity, 1LL << (max_lg_order - lg_players));
    std::vector<mpz_class> A;
    std::vector<mpz_class> eval_result;
    std::vector<bn::Ec1> B;
    int degree = 1 << lg_degree;

    for(int i = 1; i < degree; ++i)
        A.push_back(prover_coef[i]);
    for(int i = degree - 2; i >= 0; --i)
        B.push_back(pp[i]);

    std::vector<bn::Ec1> C;

    fft_eval(prover_coef, omega, lg_players, eval_result);

    //padding
    A.push_back(0);
    B.push_back(g1 * 0);

    conv(A, B, C);

    std::vector<bn::Ec1> H;
    H.resize(degree);
    for(int i = 0; i < degree - 1; ++i)
    {
        H[i] = C[i + degree - 2];
    }
    H[degree - 1] = g1 * 0;
    ret.resize(1 << lg_players);
    std::vector<bn::Ec1> pi_array;
    pi_array.resize(1 << lg_players);
    fast_fourier_transform(H.data(), degree, 1 << lg_players, omega, pi_array.data());
    for(int i = 0; i < (1 << lg_players); ++i)
    {
        ret[i].pi = pi_array[i];
        ret[i].ans = eval_result[i];
    }
    return ret;
}

#endif //DKG_VSS_KZG_IMPROVED_PROVER_H
