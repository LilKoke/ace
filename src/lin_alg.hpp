#ifndef _CLASS_LINALG_H_
#define _CLASS_LINALG_H_
#include <bits/stdc++.h>

class LinAlg
{
public:
    LinAlg();
    ~LinAlg();
    std::vector<std::vector<double>> dot(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B);
    std::vector<std::vector<double>> inv(std::vector<std::vector<double>> A);
    std::vector<std::vector<double>> T(std::vector<std::vector<double>> A);
    std::vector<std::vector<double>> T(std::vector<double> x);
    std::vector<std::vector<double>> add(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B);
    std::vector<double> add(std::vector<double> x, std::vector<double> y);
    std::vector<double> mult(double a, std::vector<double> x);
    std::vector<std::vector<double>> mult(double a, std::vector<std::vector<double>> A);
    std::vector<double> qdot(std::vector<double> p, std::vector<double> q);
    std::vector<double> squeeze(std::vector<std::vector<double>> A);
    std::vector<std::vector<double>> I(int n);
    double norm(std::vector<double> q);
};

#endif