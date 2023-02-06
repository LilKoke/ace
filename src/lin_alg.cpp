#include "lin_alg.hpp"
#include <bits/stdc++.h>

LinAlg::LinAlg(){
}

LinAlg::~LinAlg(){
}

std::vector<std::vector<double> > LinAlg::dot(std::vector<std::vector<double> > A, std::vector<std::vector<double> > B){
    std::vector<std::vector<double> > AdotB(A.size());
    for(int i = 0; i < A.size(); i++){
        for(int j = 0; j < B[0].size(); j++){
            double sum = 0;
            for(int k = 0; k < A[0].size(); k++){
                sum += A[i][k] * B[k][j];
            }
            AdotB[i].push_back(sum);
        }
    }
    return AdotB;
}

std::vector<std::vector<double> > LinAlg::inv(std::vector<std::vector<double> > A){
    std::vector<std::vector<double> > inv_A(A.size());
    for(int i=0;i<A.size();i++){
        for(int j=0;j<A.size();j++){
            inv_A[i].push_back((i==j)?1.0:0.0);
        }
    }
    for(int i=0;i<A.size();i++){
        double buf=1/A[i][i];
        assert(A[i][i] != 0);
        for(int j=0;j<A.size();j++){
            A[i][j]*=buf;
            inv_A[i][j]*=buf;
        }
        for(int j=0;j<A.size();j++){
            if(i!=j){
                buf=A[j][i];
                for(int k=0;k<A.size();k++){
                    A[j][k]-=A[i][k]*buf;
                    inv_A[j][k]-=inv_A[i][k]*buf;
                }
            }
        }
    }
    return inv_A;
}

std::vector<std::vector<double> > LinAlg::T(std::vector<std::vector<double> > A){
    std::vector<std::vector<double> > AT(A[0].size());
    for(int i = 0; i < A[0].size(); i++){
        for(int j = 0; j < A.size(); j++){
            AT[i].push_back(A[j][i]);
        }
    }
    return AT;
}

std::vector<std::vector<double> > LinAlg::add(std::vector<std::vector<double> > A, std::vector<std::vector<double> > B){
    std::vector<std::vector<double> > sumAB(A.size());
    for(int i = 0; i < A.size(); i++){
        for(int j = 0; j < A[i].size(); j++){
            sumAB[i].push_back(A[i][j] + B[i][j]);
        }
    }
    return sumAB;
}

std::vector<double> LinAlg::add(std::vector<double> x, std::vector<double> y){
    std::vector<double> sumxy(x.size());
    for(int i = 0; i < x.size(); i++){
        sumxy[i] = x[i] + y[i];
    }
    return sumxy;
}

std::vector<std::vector<double> > LinAlg::mult(double a, std::vector<std::vector<double> > A){
    std::vector<std::vector<double> > multaA(A.size());
    for(int i = 0; i < A.size(); i++){
        for(int j = 0; j < A[i].size(); j++){
            multaA[i].push_back(a * A[i][j]);
        }
    }
    return multaA;
}