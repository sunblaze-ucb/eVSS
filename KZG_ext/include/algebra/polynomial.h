//
// Created by Tiancheng Xie on 3/10/2021.
//

#ifndef DKG_VSS_POLYNOMIAL_H
#define DKG_VSS_POLYNOMIAL_H
#include "algebra/prime_field.h"
template<class coef_T, class variable_T, class eval_T>
class polynoimal
{
public:
    coef_T* data;
    int d;
    polynoimal()
    {
        data = NULL;
        d = 0;
    }
    void init(int degree, coef_T* coefficient)
    {
        d = degree;
        if(data != NULL)
            delete[] data;

        data = new T[d];
        for(int i = 0; i < degree; ++i)
            data[i] = coefficient;
    }
    eval_T eval(const variable_T &x)
    {
        eval_T ret = 0;
        variable_T pow_x = 1;
        for(int i = 0; i < d; ++i)
        {
            ret = ret + data[i] * pow_x;
            pow_x = pow_x * x;
        }
        return ret;
    }
};


#endif //DKG_VSS_POLYNOMIAL_H
