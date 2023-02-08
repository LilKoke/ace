#ifndef _CLASS_DYNAMICS_H_
#define _CLASS_DYNAMICS_H_
#include <bits/stdc++.h>
#include "lin_alg.hpp"

class Dynamics
{
    public:
    Dynamics();
    ~Dynamics();
    std::vector<double> calculate_domega(std::vector<double> omega, std::vector<double> M, std::vector<double> Iner);
    std::vector<double> calculate_dq(LinAlg la, std::vector<double> q, std::vector<double> omega);
};

#endif