#include <bits/stdc++.h>
#include "../src/lin_alg.hpp"

using namespace std;

int main(){
    cout << "test_linalg.cpp: ";
    LinAlg la;
    vector<vector<double> > x = {{1}, {2}};
    vector<vector<double> > y = {{0}, {3}};    
    vector<vector<double> > A = {{1, 2}, {3, 4}};
    vector<vector<double> > B = {{1, 2, 3}, {4, 5, 6}};
    vector<vector<double> > P = {{0, 1}, {2, 3}};
    vector<vector<double> > Q = {{0, 0.1}, {0.2, 0.3}};
    vector<vector<double> > test1, \
                            test2, \
                            test3, \
                            test4, \
                            test5, \
                            test6, \
                            test7, \
                            test8, \
                            test9, \
                            test10;
    vector<double> q = {1, 0, 0, 0};
    vector<double> omega = {0, 0.2, 0, 0};
    vector<double> dq = la.mult(0.5, la.qdot(q, omega));
    double k = 3;
    vector<vector<double> > kA = la.mult(k, A);
    vector<double> sx = la.squeeze(x);
    // dot
    test1 = la.dot(x, la.T(y));
    test2 = la.dot(la.T(x), y);
    assert(test1.size() == 2);
    assert(test1[0][1] == 3);
    assert(test1[1][1] == 6);
    cout << ".";
    assert(test2.size() == 1);
    assert(test2[0][0] == 6);
    cout << ".";
    test3 = la.dot(A, x);
    assert(test3[0][0] == 5);
    assert(test3[1][0] == 11);
    cout << ".";
    test4 = la.dot(A, P);
    assert(test4.size() == 2 && test4[0].size() == 2);
    assert(test4[0][0] == 4);
    cout << ".";
    test5 = la.dot(A, B);
    assert(test5.size() == 2 && test5[0].size() == 3);
    assert(test5[0][0] == 9);
    assert(test5[1][2] == 33);
    cout << ".";
    // inv, I
    test6 = la.inv(A);
    assert(test6[0][0] = -A[1][1]/2);
    assert(test6[0][1] = -A[1][0]/2);
    cout << ".";
    test7 = la.inv(la.I(6));
    assert(test7[0][0] == 1);
    assert(test7[3][4] == 0);
    cout << ".";
    // T
    test8 = la.T(A);
    assert(test8[0][1] == 3);
    assert(test8[1][0] == 2);
    cout << ".";
    // add mult
    test9 = la.add(A, P);
    assert(test9[0][0] == 1);
    cout << ".";
    test10 = la.add(A, la.mult(-1, Q));
    assert(test10[0][1] == 1.9);
    cout << ".";
    // qdot
    assert(dq[1] == 0.1);
    cout << ".";
    // mult
    assert(kA[0][0] == k);    
    // squeeze    
    assert(sx[0] == 1);
    assert(sx[1] == 2);
    cout << ".";
    cout << " PASSED" << endl;
}