#include "dynamics.hpp"
#include "lin_alg.hpp"
#include <bits/stdc++.h>

Dynamics::Dynamics(){
}

Dynamics::~Dynamics(){
}

std::vector<double> Dynamics::calculate_domega(std::vector<double> omega, std::vector<double> M, std::vector<double> Iner){
    std::vector<double> domega(4);
    domega[0] = 0;
    domega[1] = (M[0] - (Iner[2] - Iner[1]) * omega[2] * omega[3]) / Iner[0]; 
    domega[2] = (M[1] - (Iner[0] - Iner[2]) * omega[1] * omega[3]) / Iner[1];
    domega[3] = (M[2] - (Iner[1] - Iner[0]) * omega[1] * omega[2]) / Iner[2];
    return domega; 
}

std::vector<double> Dynamics::calculate_dq(LinAlg la, std::vector<double> q, std::vector<double> omega){
    std::vector<double> dq(4);
    dq = la.mult(0.5, la.qdot(q, omega));
    return dq;
}


