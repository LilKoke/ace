#include <bits/stdc++.h>
#include "../src/lin_alg.hpp"

using namespace std;

int main(){
    LinAlg la;
    vector<vector<double> > x = {{1}, {2}};
    vector<vector<double> > y = {{0}, {3}};
    vector<vector<double> > xxT = la.dot(x, la.T(x));
    cout << xxT[0][0] << " " << xxT[0][1] << endl;
    cout << xxT[1][0] << " " << xxT[1][1] << endl;
    vector<vector<double> > A = {{1, 2}, {3, 4}};
    vector<vector<double> > P = {{0, 1}, {2, 3}};
    vector<vector<double> > inv_A = la.inv(A);
    cout << inv_A[0][0] << " " << inv_A[0][1] << endl;
    cout << inv_A[1][0] << " " << inv_A[1][1] << endl;
    vector<vector<double> > Ax = la.dot(A, x);
    vector<vector<double> > Ay = la.dot(A, y);
    cout << Ax[0][0] << endl;
    cout << Ax[1][0] << endl;
    vector<vector<double> > sumAxAy = la.add(Ax, Ay);
    cout << sumAxAy[0][0] << endl;
    cout << sumAxAy[1][0] << endl;
    double k = 3;
    vector<vector<double> > kA = la.mult(k, A);
    cout << kA[0][0] << " " << kA[0][1] << endl;
    cout << kA[1][0] << " " << kA[1][1] << endl;
    vector<double> sx = la.squeeze(x);
    cout << sx[0] << sx[1] << endl;
    vector<double> q = {1, 0, 0, 0};
    vector<double> omega = {0, 0.2, 0, 0};
    vector<double> dq = la.mult(0.5, la.qdot(q, omega));
    assert(dq[1] == 0.1);
    double h = 0.1;
    vector<double> dq2 = la.mult(0.5*h, la.qdot(q, omega));
    cout << dq2[0] << " " << dq2[1] << " " << dq2[2] << " " << dq2[3] << endl;
    return 0;
}