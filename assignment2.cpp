#include <bits/stdc++.h>
#include <sys/stat.h>
#include "src/lin_alg.hpp"

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
    double omega_nominal = 17 * 2 * PI / 60;
    double omega[4] = {0, 0.1, omega_nominal + 0.1, 0};
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
    ofstream writing_file1, writing_file2, writing_file3, writing_file4;
    string folder_name = "out/assignment2";
    mkdir(folder_name.c_str(), 0777);
    string filename_omega = "out/assignment2/omega.csv";
    string filename_q = "out/assignment2/q.csv";
    string filename_domega = "out/assignment2/domega.csv";
    string filename_dq = "out/assignment2/dq.csv";
    writing_file1.open(filename_omega, ios::out);
    writing_file2.open(filename_q, ios::out);
    writing_file3.open(filename_domega, ios::out);
    writing_file4.open(filename_dq, ios::out);
    writing_file1 << "t" << " " << "omega1" << " " << "omega2" << " " << "omega3" << endl;
    writing_file2 << "t" << " " << "q0" << " " << "q1" << " " << "q2" << " " << "q3" << endl;
    writing_file3 << "t" << " " << "domega1" << " " << "domega2" << " " << "domega3" << endl;
    writing_file4 << "t" << " " << "dq0" << " " << "dq1" << " " << "dq2" << " " << "dq3" << endl;
    vector<vector<double> > x = {{0}, {0}, {0}, {0}, {0}, {0}, {0}};
    vector<vector<double> > z;
    vector<vector<double> > x_bar, x_hat; x_hat = x;
    vector<vector<double> > A, B, H, K, covM, P;
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    double mean_w = 0;
    double mean_v = 0;
    double sigma_w = 0.01;
    double sigma_v = 0.01;
    double mean_s = 0;
    double sigma_s = 1;
    normal_distribution<> dist_w(mean_w, sigma_w), dist_v(mean_v, sigma_v), dist_s(mean_s, sigma_s);
    vector<double> q_s(4), omega_s(4);
    for(int i = 0; i < 4; i++) q_s[i] = dist_s(engine);
    omega_s[0] = 0;
    for(int i = 0; i < 3; i++) omega_s[i+1] = dist_s(engine);

    vector<vector<double> > w(3), v(3);
    vector<vector<double> > R = {
                                    {sigma_v, 0, 0},
                                    {0, sigma_v, 0},
                                    {0, 0, sigma_v}
                                };
    LinAlg la;
    
    while(t <= tn){
        // propagation phase
        A = {
                {0, -omega[1]/2, -omega[2]/2, -omega[3]/2, -q[1]/2, -q[2]/2, -q[3]/2},
                {omega[1]/2, 0, omega[3]/2, -omega[2]/2, q[0]/2, -q[3]/2, q[2]/2},
                {omega[2]/2, -omega[3]/2, 0, omega[1]/2, q[3]/2, q[0]/2, -q[1]/2},
                {omega[3]/2, omega[2]/2, -omega[1]/2, 0, -q[2]/2, q[1]/2, q[0]/2},
                {0, 0, 0, 0, 0, (Iner[1]-Iner[2])/Iner[0]*omega[3], (Iner[1]-Iner[2])/Iner[0]*omega[2]},
                {0, 0, 0, 0, (Iner[2]-Iner[0])/Iner[1]*omega[3], 0, (Iner[2]-Iner[0])/Iner[1]*omega[1]},
                {0, 0, 0, 0, (Iner[0]-Iner[1])/Iner[2]*omega[2], (Iner[0]-Iner[1])/Iner[2]*omega[1], 0}
            };
        B = {
                {0, 0, 0}, 
                {0, 0, 0}, 
                {0, 0, 0}, 
                {0, 0, 0}, 
                {1/Iner[0], 0, 0}, 
                {0, 1/Iner[1], 0}, 
                {0, 0, 1/Iner[2]}
            };
        P = {
                {0.1, 0, 0, 0, 0, 0, 0},
                {0, 0.1, 0, 0, 0, 0, 0},
                {0, 0, 0.1, 0, 0, 0, 0},
                {0, 0, 0, 0.1, 0, 0, 0},
                {0, 0, 0, 0, 0.1, 0, 0},
                {0, 0, 0, 0, 0, 0.1, 0},
                {0, 0, 0, 0, 0, 0, 0.1}
            };
        for(int i = 0; i < 3; i++){
            w[i].push_back(dist_w(engine));
            v[i].push_back(dist_v(engine));
            M[i] = w[i][0];
        }
        x_bar = la.add(la.dot(A, x_hat), la.dot(B, w));
        double varv = sigma_v * sigma_v;
        covM = la.add(la.dot(la.dot(A, P), la.T(A)), la.mult(varv, la.dot(B, la.T(B))));
        
        // runge-kutta
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

        // measurement
        H = {
                {2*q[0], 2*q[1], -2*q[2], -2*q[3], 0, 0, 0},
                {2*q[3], 2*q[2], 2*q[1], 2*q[0], 0, 0, 0},
                {-2*q[0], 2*q[3], -2*q[0], 2*q[1], 0, 0, 0},   
            };
        vector<double> dq(4), domega(4), dq_hat(4), domega_hat(4);
        for(int i = 0; i < 4; i++){
            dq[i] = q[i] - q_s[i];
            domega[i] = omega[i] - omega_s[i];
        }
        for(int i = 0; i < 7; i++){
            if(i <= 3) x[i][0] = dq[i];
            else x[i][0] = domega[i - 3];
        }
        z = la.add(la.dot(H, x), v);


        // calculate Kalman gain
        K = la.dot(la.dot(P, la.T(H)), la.inv(R));

        vector<vector<double> > error = la.add(z, la.mult(-1, la.dot(H, x_bar)));

        // update x_hat and P
        x_hat = la.add(x_bar, la.dot(K, error));
        P = la.inv(la.add(la.inv(covM), la.dot(la.dot(la.T(H), la.inv(R)), H)));

        for(int i = 0; i < 7; i++){
            if(i <= 3) q_s[i] += x_hat[i][0];
            else omega_s[i - 3] += x_hat[i][0];
            x[i][0] = 0; 
        }

        writing_file1 << t << " " << omega[1] << " " << omega[2] << " " << omega[3] << endl;        
        writing_file2 << t << " " << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << endl;
        writing_file3 << t << " " << omega[1] - omega_s[1] << " " << omega[2] - omega_s[2] << " " << omega[3] - omega_s[3] << endl;
        writing_file4 << t << " " << q[0] - q_s[0] << " " << q[1] - q_s[1] << " " << q[2] - q_s[2] << " " << q[3] - q_s[3] << endl;
    }
    writing_file1.close();
    writing_file2.close();
    writing_file3.close();
    writing_file4.close();
    return 0;
}