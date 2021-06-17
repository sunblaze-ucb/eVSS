//
// Created by Tiancheng Xie on 3/10/2021.
//

#ifndef DKG_VSS_ENVIROMENT_H
#define DKG_VSS_ENVIROMENT_H

#include "bn.h"
#include "test_point.hpp"

bn::Ec1 g1;
bn::Ec2 g2;
mpz_class p;
gmp_randstate_t r_state;

mpz_class root_of_unity;
const int max_lg_order = 28;
class single_proof
{
public:
    mpz_class ans;
    bn::Ec1 pi;
};
template<class exp_T>
mpz_class fastpow(mpz_class x, exp_T y)
{
    mpz_class ret = 1;
    exp_T two = 2;
    while(y != 0)
    {
        if(y % 2 == 1)
        {
            ret = ret * x % p;
        }
        x = x * x % p;
        y = y / two;
    }
    return ret;
}

void init_enviroment()
{
    p.set_str("21888242871839275222246405745257275088548364400416034343698204186575808495617",10);
    root_of_unity.set_str("12237580043724246686314448501040597564946944530769647865992922846663565257669", 10);
    auto seed = rand();
    gmp_randinit_default(r_state);
    gmp_randseed_ui(r_state, seed);


    //bilinear g1 g2
    bn::CurveParam cp = bn::CurveSNARK1;
    bn::Param::init(cp);
    const Point& pt = selectPoint(cp);
    g2 = bn::Ec2(
            bn::Fp2(bn::Fp(pt.g2.aa), bn::Fp(pt.g2.ab)),
            bn::Fp2(bn::Fp(pt.g2.ba), bn::Fp(pt.g2.bb))
    );
    g1 = bn::Ec1(pt.g1.a, pt.g1.b);
}
#endif //DKG_VSS_ENVIROMENT_H
