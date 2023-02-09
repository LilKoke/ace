#include "lin_alg.hpp"
#include <bits/stdc++.h>

LinAlg::LinAlg()
{
}

LinAlg::~LinAlg()
{
}

std::vector<std::vector<double>> LinAlg::dot(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B)
{
    std::vector<std::vector<double>> AdotB(A.size());
    for (int i = 0; i < A.size(); i++)
    {
        for (int j = 0; j < B[0].size(); j++)
        {
            double sum = 0;
            for (int k = 0; k < A[0].size(); k++)
            {
                sum += A[i][k] * B[k][j];
            }
            AdotB[i].push_back(sum);
        }
    }
    return AdotB;
}

std::vector<std::vector<double>> LinAlg::inv(std::vector<std::vector<double>> A)
{
    std::vector<std::vector<double>> inv_A(A.size());
    for (int i = 0; i < A.size(); i++)
    {
        for (int j = 0; j < A.size(); j++)
        {
            inv_A[i].push_back((i == j) ? 1.0 : 0.0);
        }
    }
    for (int i = 0; i < A.size(); i++)
    {
        double buf = 1 / A[i][i];
        assert(A[i][i] != 0);
        for (int j = 0; j < A.size(); j++)
        {
            A[i][j] *= buf;
            inv_A[i][j] *= buf;
        }
        for (int j = 0; j < A.size(); j++)
        {
            if (i != j)
            {
                buf = A[j][i];
                for (int k = 0; k < A.size(); k++)
                {
                    A[j][k] -= A[i][k] * buf;
                    inv_A[j][k] -= inv_A[i][k] * buf;
                }
            }
        }
    }
    return inv_A;
}

std::vector<std::vector<double>> LinAlg::T(std::vector<std::vector<double>> A)
{
    std::vector<std::vector<double>> AT(A[0].size());
    for (int i = 0; i < A[0].size(); i++)
    {
        for (int j = 0; j < A.size(); j++)
        {
            AT[i].push_back(A[j][i]);
        }
    }
    return AT;
}

std::vector<std::vector<double>> LinAlg::T(std::vector<double> x)
{
    std::vector<std::vector<double>> xT(x.size());
    for (int i = 0; i < x.size(); i++)
    {
        xT[i].push_back(x[i]);
    }
    return xT;
}

std::vector<std::vector<double>> LinAlg::add(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B)
{
    std::vector<std::vector<double>> sumAB(A.size());
    for (int i = 0; i < A.size(); i++)
    {
        for (int j = 0; j < A[i].size(); j++)
        {
            sumAB[i].push_back(A[i][j] + B[i][j]);
        }
    }
    return sumAB;
}

std::vector<double> LinAlg::add(std::vector<double> x, std::vector<double> y)
{
    std::vector<double> sumxy(x.size());
    for (int i = 0; i < x.size(); i++)
    {
        sumxy[i] = x[i] + y[i];
    }
    return sumxy;
}

std::vector<double> LinAlg::mult(double a, std::vector<double> x)
{
    std::vector<double> multax(x.size());
    for (int i = 0; i < x.size(); i++)
        multax[i] = x[i] * a;
    return multax;
}

std::vector<std::vector<double>> LinAlg::mult(double a, std::vector<std::vector<double>> A)
{
    std::vector<std::vector<double>> multaA(A.size());
    for (int i = 0; i < A.size(); i++)
    {
        for (int j = 0; j < A[i].size(); j++)
        {
            multaA[i].push_back(a * A[i][j]);
        }
    }
    return multaA;
}

std::vector<double> LinAlg::qdot(std::vector<double> p, std::vector<double> q)
{
    std::vector<double> pqdotq;
    assert(p.size() == 4);
    assert(q.size() == 4);
    pqdotq = {
        {p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3]},
        {p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2]},
        {p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1]},
        {p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0]}};
    return pqdotq;
}

std::vector<double> LinAlg::squeeze(std::vector<std::vector<double>> A)
{
    std::vector<double> x(A.size());
    assert(A[0].size() == 1);
    for (int i = 0; i < A.size(); i++)
        x[i] = A[i][0];
    return x;
}

std::vector<std::vector<double>> LinAlg::I(int n)
{
    std::vector<std::vector<double>> Im(n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
                Im[i].push_back(1.0);
            else
                Im[i].push_back(0.0);
        }
    }
    return Im;
}

double LinAlg::norm(std::vector<double> q)
{
    double sum = 0;
    for (int i = 0; i < q.size(); i++)
    {
        sum += q[i] * q[i];
    }
    return sqrt(sum);
}