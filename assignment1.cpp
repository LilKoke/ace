#include <bits/stdc++.h>
#include <sys/stat.h>
#include "src/lin_alg.hpp"

using namespace std;

LinAlg la;

vector<double> calculate_domega(vector<double> omega, vector<double> M, vector<double> Iner){
    vector<double> domega(4);
    domega[0] = 0;
    domega[1] = (M[0] - (Iner[2] - Iner[1]) * omega[2] * omega[3]) / Iner[0]; 
    domega[2] = (M[1] - (Iner[0] - Iner[2]) * omega[1] * omega[3]) / Iner[1];
    domega[3] = (M[2] - (Iner[1] - Iner[0]) * omega[1] * omega[2]) / Iner[2];
    return domega; 
}

// vector<double> calculate_dq(vector<double> q, vector<double> omega){
//     vector<double> dq(4);
//     LinAlg la;
//     vector<vector<double> > Omega = {
//         {0, -omega[1], -omega[2], -omega[3]},
//         {omega[1], 0, omega[3], -omega[2]},
//         {omega[2], -omega[3], 0, omega[1]},
//         {omega[3], omega[2], -omega[1], 0}
//     };
//     dq = la.squeeze(la.dot(Omega, la.T(q)));
//     return dq;
// }

vector<double> calculate_dq(vector<double> q, vector<double> omega){
    vector<double> dq(4);
    dq = la.mult(0.5, la.qdot(q, omega));
    return dq;
}

// vector<double> calculate_dq(vector<double> q, vector<double> omega){
//     vector<double> dq(4);
//     vector<vector<double> > qmat = {
//         {0, -q[1], -q[2], -q[3]},
//         {0, q[0], -q[3], q[2]},
//         {0, q[3], q[0], -q[1]},
//         {0, -q[2], q[1], q[0]}
//     };
//     dq = la.squeeze(la.mult(0.5, la.dot(qmat, la.T(omega))));
//     return dq;
// }


#define PI 3.1415

int main(){
    vector<double> q = {1, 0, 0, 0};
    vector<double> dq;
    double omega_s = 17 * 2 * PI / 60;
    vector<double> omega = {0, 0.1, omega_s + 0.1, 0};
    vector<double> M = {0, 0, 0};
    vector<double> Iner = {1.9, 1.6, 2.0};
    vector<double> domega;
    double t0 = 0;
    double tn = 30;
    double n = tn * 1000;
    double h = (tn - t0)/n;
    vector<double> k1_omega(4), k2_omega(4), k3_omega(4), k4_omega(4);
    vector<double> k1_q(4), k2_q(4), k3_q(4), k4_q(4);
    vector<double> omega_runge(4), q_runge(4);
    double t = 0;
    ofstream writing_file1, writing_file2;
    string folder_name = "out/assignment1";
    mkdir(folder_name.c_str(), 0777);
    string filename_omega = "out/assignment1/omega.csv";
    string filename_q = "out/assignment1/q.csv";
    writing_file1.open(filename_omega, ios::out);
    writing_file2.open(filename_q, ios::out);
    writing_file1 << "t" << " " << "omega1" << " " << "omega2" << " " << "omega3" << endl;
    writing_file2 << "t" << " " << "q0" << " " << "q1" << " " << "q2" << " " << "q3" << endl;
    while(t <= tn){
        for(int i = 0; i < 4; i++){
            omega_runge[i] = omega[i];
            q_runge[i] = q[i];
        }
        domega = calculate_domega(omega_runge, M, Iner);
        dq = calculate_dq(q_runge, omega_runge);
        k1_omega = la.mult(h, domega);
        k1_q = la.mult(h, dq);
        omega_runge = la.add(omega_runge, la.mult(0.5, k1_omega));
        q_runge = la.add(q_runge, la.mult(0.5, k1_q));

        domega = calculate_domega(omega_runge, M, Iner);
        dq = calculate_dq(q_runge, omega_runge);
        k2_omega = la.mult(h, domega);
        k2_q = la.mult(h, dq);
        omega_runge = la.add(omega_runge, la.mult(0.5, k2_omega));
        q_runge = la.add(q_runge, la.mult(0.5, k2_q));

        domega = calculate_domega(omega_runge, M, Iner);
        dq = calculate_dq(q_runge, omega_runge);
        k3_omega = la.mult(h, domega);
        k3_q = la.mult(h, dq);
        omega_runge = la.add(omega_runge, k3_omega);
        q_runge = la.add(q_runge, la.mult(1, k3_q));

        domega = calculate_domega(omega_runge, M, Iner);
        dq = calculate_dq(q_runge, omega_runge);
        k4_omega = la.mult(h, domega);
        k4_q = la.mult(h, dq);

        for(int i = 0; i < 4; i++){
            omega[i] += (k1_omega[i] + 2 * k2_omega[i] + 2 * k3_omega[i] + k4_omega[i])/6;
            q[i] += (k1_q[i] + (2 * k2_q[i]) + (2 * k3_q[i]) + k4_q[i])/6;
        }
        t += h;
        writing_file1 << t << " " << omega[1] << " " << omega[2] << " " << omega[3] << endl;        
        writing_file2 << t << " " << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << endl;
        
    }
    writing_file1.close();
    writing_file2.close();
    return 0;
}