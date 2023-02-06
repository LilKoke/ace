#ifndef _CLASS_LINALG_H_
#define _CLASS_LINALG_H_
#include <bits/stdc++.h>

class LinAlg
{
    public:
    LinAlg();
    ~LinAlg();
    std::vector<std::vector<double> > dot(std::vector<std::vector<double> > A, std::vector<std::vector<double> > B);
    std::vector<std::vector<double> > inv(std::vector<std::vector<double> > A);
    std::vector<std::vector<double> > T(std::vector<std::vector<double> > A);
    std::vector<std::vector<double> > add(std::vector<std::vector<double> > A, std::vector<std::vector<double> > B);
    std::vector<double> add(std::vector<double> x, std::vector<double> y);
    std::vector<std::vector<double> > mult(double, std::vector<std::vector<double> > A);
};

#endif