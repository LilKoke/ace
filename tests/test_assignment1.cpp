#include <bits/stdc++.h>

using namespace std;

void eom_rotate(double omega[], double M[], double Iner[], double domega[]){
    domega[0] = (Iner[2] - Iner[1]) * omega[1] * omega[2] - M[0]; 
    domega[1] = (Iner[0] - Iner[2]) * omega[0] * omega[2] - M[0];
    domega[2] = (Iner[1] - Iner[0]) * omega[0] * omega[1] - M[0]; 
}

void calculate_dq(double q[], double omega[], double dq[]){
    dq[0] = (q[0] * omega[0] - q[1] * omega[1] - q[2] * omega[2] - q[3] * omega[3]) / 2;
    dq[1] = (q[0] * omega[1] + q[1] * omega[0] + q[2] * omega[3] - q[3] * omega[2]) / 2;
    dq[2] = (q[0] * omega[2] - q[1] * omega[3] + q[2] * omega[0] + q[3] * omega[1]) / 2;
    dq[3] = (q[0] * omega[3] + q[1] * omega[2] - q[2] * omega[1] + q[3] * omega[0]) / 2;
}

void update_q(double q[], double dq[], double dt){
    double q_sum = 0;
    for(int i = 0; i < 4; i++){
        q[i] += dq[i] * dt;
    }
    for(int i = 0; i < 4; i++) q_sum += q[i];
    for(int i = 0; i < 4; i++) q[i] /= q_sum;
}
#define PI 3.1415

int main(){
    double q[4] = {1, 0, 0, 0};
    double dq[4] = {0, 0, 0, 0};
    double omega[4] = {0, 0.2, 0, 0};
    double domega[3];
    double dt = 1;

    calculate_dq(q, omega, dq);
    update_q(q, dq, dt);
    assert(q[0] + q[1] + q[2] + q[3] - 1 < 1e-10);
    assert(q[0] == 1 / 1.1);
    assert(q[1] == 0.1 / 1.1);
    cout << "unit_test.cpp: " << "PASSED" << endl;
    return 0;
}