#include <bits/stdc++.h>

using namespace std;

void calculate_domega(double omega[], double M[], double Iner[], double domega[]){
    domega[0] = 0;
    domega[1] = ((Iner[2] - Iner[1]) * omega[1] * omega[2] - M[0]) / Iner[0]; 
    domega[2] = ((Iner[0] - Iner[2]) * omega[0] * omega[2] - M[0]) / Iner[1];
    domega[3] = ((Iner[1] - Iner[0]) * omega[0] * omega[1] - M[0]) / Iner[2]; 
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
    double omega_s = 17 * 2 * PI / 60;
    double omega[4] = {0, 0.1, omega_s + 0.1, 0};
    double M[3] = {0, 0, 0};
    double Iner[3] = {1.9, 1.6, 2.0};
    double domega[3];
    double t0 = 0;
    double tn = 1;
    double n = 100;
    double h = (tn - t0) / n;
    double k1_omega[4], k2_omega[4], k3_omega[4], k4_omega[4];
    double k1_q[4], k2_q[4], k3_q[4], k4_q[4];
    double omega_runge[4], q_runge[4];
    double t = 0;
    ofstream writing_file1, writing_file2;
    string filename_omega = "out/omega.csv";
    string filename_q = "out/q.csv";
    writing_file1.open(filename_omega, ios::out);
    writing_file2.open(filename_q, ios::out);
    writing_file1 << "t" << " " << "omega1" << " " << "omega2" << " " << "omega3" << endl;
    writing_file2 << "t" << " " << "q0" << " " << "q1" << " " << "q2" << " " << "q3" << endl;
    while(t <= tn){
        for(int i = 0; i < 4; i++){
            omega_runge[i] = omega[i];
            q_runge[i] = q[i];
        }
        calculate_domega(omega_runge, M, Iner, domega);
        calculate_dq(q_runge, omega_runge, dq);
        for(int i = 0; i < 4; i++){
            k1_omega[i] = h * domega[i];
            k1_q[i] = h * dq[i];
            omega_runge[i] += k1_omega[i] / 2;
            q_runge[i] += k1_q[i] / 2;
        }
        calculate_domega(omega_runge, M, Iner, domega);
        calculate_dq(q_runge, omega_runge, dq);
        for(int i = 0; i < 4; i++){
            k2_omega[i] = h * domega[i];
            k2_q[i] = h * dq[i];
            omega_runge[i] += k2_omega[i] / 2;
            q_runge[i] += k2_q[i] / 2;
        }
        calculate_domega(omega_runge, M, Iner, domega);
        calculate_dq(q_runge, omega_runge, dq);
        for(int i = 0; i < 4; i++){
            k3_omega[i] = h * domega[i];
            k3_q[i] = h * dq[i];
            omega_runge[i] += k3_omega[i];
            q_runge[i] += k3_q[i];
        }
        calculate_domega(omega_runge, M, Iner, domega);
        calculate_dq(q_runge, omega_runge, dq);
        for(int i = 0; i < 4; i++){
            k4_omega[i] = h * domega[i];
            k4_q[i] = h * dq[i];
        }

        for(int i = 0; i < 4; i++){
            omega[i] += k1_omega[i] + 2 * k2_omega[i] + 2 * k3_omega[i] + k4_omega[i];
            q[i] += k1_q[i] + 2 * k2_q[i] + 2 * k3_q[i] + k4_q[i];
        }
        t += h;
        writing_file1 << t << " " << omega[1] << " " << omega[2] << " " << omega[3] << endl;        
        writing_file2 << t << " " << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << endl;
        
    }
    writing_file1.close();
    writing_file2.close();
    return 0;
}