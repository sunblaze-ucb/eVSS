//
// Created by Tiancheng Xie on 3/10/2021.
//

#ifndef DKG_VSS_KZG_VERIFIER_H
#define DKG_VSS_KZG_VERIFIER_H

#include "bn.h"
#include "test_point.hpp"
#include "gmpxx.h"
#include "gmp.h"
//verify \phi(x) == y
bool verify(const mpz_class &x, const mpz_class &y, const bn::Ec1 &pi, const bn::Ec1 &c, const bn::Ec1 g1, const bn::Ec2 g2, const bn::Ec2 g_tau)
{
    bn::Fp12 e1, e2;
    auto e1x = g2;
    bn::Ec1 e1y = c - (g1 * y);
    bn::Ec2 e2x = g_tau - (g2 * x);
    auto e2y = pi;

    bn::opt_atePairing(e1, e1x, e1y);
    bn::opt_atePairing(e2, e2x, e2y);
    return e1 == e2;
}

#endif //DKG_VSS_KZG_VERIFIER_H
