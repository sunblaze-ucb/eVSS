//
// Created by Tiancheng Xie on 3/12/2021.
//

#ifndef DKG_VSS_KZG_COMMIT_H
#define DKG_VSS_KZG_COMMIT_H

std::vector<mpz_class> prover_coef;
bn::Ec1 commit_proof(int lg_degree, const std::vector<bn::Ec1> &pp)
{
    int degree = 1 << lg_degree;
    prover_coef.resize(degree);
    for(int i = 0; i < degree; ++i)
    {
        mpz_urandomm(prover_coef[i].get_mpz_t(), r_state, p.get_mpz_t());
    }
    bn::Ec1 ret = g1 * 0;
    for(int i = 0; i < degree; ++i)
    {
        ret = ret + pp[i] * prover_coef[i];
    }
    printf("commit done\n");
    return ret;
}

#endif //DKG_VSS_KZG_COMMIT_H
