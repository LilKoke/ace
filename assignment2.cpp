#include <bits/stdc++.h>
#include <sys/stat.h>
#include "src/lin_alg.hpp"
#include "src/dynamics.hpp"

using namespace std;

LinAlg la;
Dynamics d;

#define PI 3.1415

int main(){
    vector<double> q = {1, 0, 0, 0};
    vector<double> dq;
    double omega_nominal = 17 * 2 * PI / 60;
    vector<double> omega = {0, 0.1, omega_nominal + 0.1, 0};
    vector<double> M = {0, 0, 0};
    vector<double> M_s = {0, 0, 0};
    vector<double> Iner = {1.9, 1.6, 2.0};
    vector<double> domega;
    double t0 = 0;
    double tn = 30;
    double dt = 0.01;
    vector<double> k1_omega(4), k2_omega(4), k3_omega(4), k4_omega(4);
    vector<double> k1_q(4), k2_q(4), k3_q(4), k4_q(4);
    vector<double> omega_runge(4), q_runge(4);
    double t = 0;
    ofstream writing_file1, writing_file2, writing_file3, writing_file4, writing_file5, writing_file6;
    string folder_name = "out/assignment2";
    mkdir(folder_name.c_str(), 0777);
    string filename_omega = "out/assignment2/omega.csv";
    string filename_q = "out/assignment2/q.csv";
    string filename_domega = "out/assignment2/domega.csv";
    string filename_dq = "out/assignment2/dq.csv";
    string filename_omega_s = "out/assignment2/omega_s.csv";
    string filename_q_s = "out/assignment2/q_s.csv";
    writing_file1.open(filename_omega, ios::out);
    writing_file2.open(filename_q, ios::out);
    writing_file3.open(filename_domega, ios::out);
    writing_file4.open(filename_dq, ios::out);
    writing_file5.open(filename_omega_s, ios::out);
    writing_file6.open(filename_q_s, ios::out);
    writing_file1 << "t" << " " << "omega1" << " " << "omega2" << " " << "omega3" << endl;
    writing_file2 << "t" << " " << "q0" << " " << "q1" << " " << "q2" << " " << "q3" << endl;
    writing_file3 << "t" << " " << "domega1" << " " << "domega2" << " " << "domega3" << endl;
    writing_file4 << "t" << " " << "dq0" << " " << "dq1" << " " << "dq2" << " " << "dq3" << endl;
    writing_file5 << "t" << " " << "omega1" << " " << "omega2" << " " << "omega3" << endl;
    writing_file6 << "t" << " " << "q0" << " " << "q1" << " " << "q2" << " " << "q3" << endl;
    vector<vector<double> > x = {{0}, {0}, {0}, {0}, {0}, {0}, {0}};
    vector<vector<double> > z;
    vector<vector<double> > x_bar, x_hat; x_hat = x;
    vector<vector<vector<double> > > H_list, b_list, b_slist;
    vector<vector<double> > A, B, H, K, covM, Phi, Gamma;
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    double mean_w = 0;
    double mean_v = 0;
    double sigma_w = 0.01;
    double sigma_v = 0.01;
    double mean_s = 0;
    double sigma_s = 0.1;
    normal_distribution<> dist_w(mean_w, sigma_w), dist_v(mean_v, sigma_v), dist_s(mean_s, sigma_s);
    vector<double> q_s(4), omega_s(4);
    omega_s[0] = 0;
    for(int i = 0; i < 3; i++) omega_s[i+1] = omega[i+1] + dist_s(engine);
    for(int i = 0; i < 4; i++) q_s[i] = q[i] + dist_s(engine);

    vector<vector<double> > w(3), v(3), b, b_s;
    vector<vector<double> > Q = {
                                    {sigma_w * sigma_w, 0, 0},
                                    {0, sigma_w * sigma_w, 0},
                                    {0, 0, sigma_w * sigma_w}
                                };
    vector<vector<double> > R = {
                                    {sigma_v * sigma_v, 0, 0},
                                    {0, sigma_v * sigma_v, 0},
                                    {0, 0, sigma_v * sigma_v}
                                };
    vector<vector<double> > P = la.mult(1, la.I(7));
    int idx = 0;
    double kf_period_s = 1;
    int kf_period_cycle = int(kf_period_s / dt);
    t = 0;
    
    while(t <= tn){
        writing_file1 << t << " " << omega[1] << " " << omega[2] << " " << omega[3] << endl;        
        writing_file2 << t << " " << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << endl;
        writing_file3 << t << " " << omega[1] - omega_s[1] << " " << omega[2] - omega_s[2] << " " << omega[3] - omega_s[3] << endl;
        writing_file4 << t << " " << q[0] - q_s[0] << " " << q[1] - q_s[1] << " " << q[2] - q_s[2] << " " << q[3] - q_s[3] << endl;
        writing_file5 << t << " " << omega_s[1] << " " << omega_s[2] << " " << omega_s[3] << endl;        
        writing_file6 << t << " " << q_s[0] << " " << q_s[1] << " " << q_s[2] << " " << q_s[3] << endl;


        // calculate true value of omega and q
        omega_runge = omega;
        q_runge = q;
        for(int i = 0; i < 3; i++){
            w[i].push_back(dist_w(engine));
            v[i].push_back(dist_v(engine));
            M[i] = w[i][0];
        }
        domega = d.calculate_domega(omega_runge, M, Iner);
        dq = d.calculate_dq(la, q_runge, omega_runge);
        k1_omega = la.mult(dt, domega);
        k1_q = la.mult(dt, dq);
        omega_runge = la.add(omega_runge, la.mult(0.5, k1_omega));
        q_runge = la.add(q_runge, la.mult(0.5, k1_q));

        domega = d.calculate_domega(omega_runge, M, Iner);
        dq = d.calculate_dq(la, q_runge, omega_runge);
        k2_omega = la.mult(dt, domega);
        k2_q = la.mult(dt, dq);
        omega_runge = la.add(omega_runge, la.mult(0.5, k2_omega));
        q_runge = la.add(q_runge, la.mult(0.5, k2_q));

        domega = d.calculate_domega(omega_runge, M, Iner);
        dq = d.calculate_dq(la, q_runge, omega_runge);
        k3_omega = la.mult(dt, domega);
        k3_q = la.mult(dt, dq);
        omega_runge = la.add(omega_runge, k3_omega);
        q_runge = la.add(q_runge, la.mult(1, k3_q));

        domega = d.calculate_domega(omega_runge, M, Iner);
        dq = d.calculate_dq(la, q_runge, omega_runge);
        k4_omega = la.mult(dt, domega);
        k4_q = la.mult(dt, dq);

        for(int i = 0; i < 4; i++){
            omega[i] += (k1_omega[i] + 2 * k2_omega[i] + 2 * k3_omega[i] + k4_omega[i])/6;
            q[i] += (k1_q[i] + (2 * k2_q[i]) + (2 * k3_q[i]) + k4_q[i])/6;
        }
        q = d.normalize_q(la, q);

        // simulate omega_s and q_s
        omega_runge = omega_s;
        q_runge = q_s;
        domega = d.calculate_domega(omega_runge, M_s, Iner);
        dq = d.calculate_dq(la, q_runge, omega_runge);
        k1_omega = la.mult(dt, domega);
        k1_q = la.mult(dt, dq);
        omega_runge = la.add(omega_runge, la.mult(0.5, k1_omega));
        q_runge = la.add(q_runge, la.mult(0.5, k1_q));

        domega = d.calculate_domega(omega_runge, M_s, Iner);
        dq = d.calculate_dq(la, q_runge, omega_runge);
        k2_omega = la.mult(dt, domega);
        k2_q = la.mult(dt, dq);
        omega_runge = la.add(omega_runge, la.mult(0.5, k2_omega));
        q_runge = la.add(q_runge, la.mult(0.5, k2_q));

        domega = d.calculate_domega(omega_runge, M_s, Iner);
        dq = d.calculate_dq(la, q_runge, omega_runge);
        k3_omega = la.mult(dt, domega);
        k3_q = la.mult(dt, dq);
        omega_runge = la.add(omega_runge, k3_omega);
        q_runge = la.add(q_runge, k3_q);

        domega = d.calculate_domega(omega_runge, M_s, Iner);
        dq = d.calculate_dq(la, q_runge, omega_runge);
        k4_omega = la.mult(dt, domega);
        k4_q = la.mult(dt, dq);

        for(int i = 0; i < 4; i++){
            omega_s[i] += (k1_omega[i] + 2 * k2_omega[i] + 2 * k3_omega[i] + k4_omega[i])/6;
            q_s[i] += (k1_q[i] + (2 * k2_q[i]) + (2 * k3_q[i]) + k4_q[i])/6;
        }
        q_s = d.normalize_q(la, q_s);

        A = {
                {0, -omega_s[1]/2, -omega_s[2]/2, -omega_s[3]/2, -q_s[1]/2, -q_s[2]/2, -q_s[3]/2},
                {omega_s[1]/2, 0, omega_s[3]/2, -omega_s[2]/2, q_s[0]/2, -q_s[3]/2, q_s[2]/2},
                {omega_s[2]/2, -omega_s[3]/2, 0, omega_s[1]/2, q_s[3]/2, q_s[0]/2, -q_s[1]/2},
                {omega_s[3]/2, omega_s[2]/2, -omega_s[1]/2, 0, -q_s[2]/2, q_s[1]/2, q_s[0]/2},
                {0, 0, 0, 0, 0, (Iner[1]-Iner[2])/Iner[0]*omega_s[3], (Iner[1]-Iner[2])/Iner[0]*omega_s[2]},
                {0, 0, 0, 0, (Iner[2]-Iner[0])/Iner[1]*omega_s[3], 0, (Iner[2]-Iner[0])/Iner[1]*omega_s[1]},
                {0, 0, 0, 0, (Iner[0]-Iner[1])/Iner[2]*omega_s[2], (Iner[0]-Iner[1])/Iner[2]*omega_s[1], 0}
            };
        Phi = la.add(la.I(A.size()), la.mult(dt, A));
        B = {
                {0, 0, 0}, 
                {0, 0, 0}, 
                {0, 0, 0}, 
                {0, 0, 0}, 
                {1/Iner[0], 0, 0}, 
                {0, 1/Iner[1], 0}, 
                {0, 0, 1/Iner[2]}
            };
        Gamma = la.mult(dt, B);
        covM = la.add(la.dot(Phi, la.dot(P, la.T(Phi))), la.dot(Gamma, la.dot(Q, la.T(Gamma))));
        assert(la.dot(Phi, la.dot(P, la.T(Phi))).size() == la.dot(Gamma, la.dot(Q, la.T(Gamma))).size());
        P = covM;

        // equation of state
        if(idx % kf_period_cycle == 0){
            
            // propagation phase
            x_bar = la.add(la.dot(Phi, x_hat), la.dot(Gamma, w));
            assert(la.dot(Phi, x_hat).size() == la.dot(Gamma, w).size());

            for(int i = 0; i < 3; i++){
                w[i][0] = dist_w(engine);
                v[i][0] = dist_v(engine);
            }
            
            // measurement
            H_list = {
                    {
                        {2*q_s[0], 2*q_s[1], -2*q_s[2], -2*q_s[3], 0, 0, 0},
                        {2*q_s[3], 2*q_s[2], 2*q_s[1], 2*q_s[0], 0, 0, 0},
                        {-2*q_s[0], 2*q_s[3], -2*q_s[0], 2*q_s[1], 0, 0, 0}
                    },
                };
            H = H_list[0];
            for(int i = 0; i < 7; i++){
                if(i <= 3){
                    x[i][0] = q[i] - q_s[i];
                }
                else{
                    x[i][0] = omega[i - 3] - omega_s[i - 3];
                }
            }
            // b_list = {
            //     {
            //         {q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3]},
            //         {2*(q[1]*q[2] + q[0]*q[3])},
            //         {2*(q[1]*q[3] - q[0]*q[2])}
            //     }
            // };
            
            // b_slist = {
            //     {
            //         {q_s[0]*q_s[0] + q_s[1]*q_s[1] - q_s[2]*q_s[2] - q_s[3]*q_s[3]},
            //         {2*(q_s[1]*q_s[2] + q_s[0]*q_s[3])},
            //         {2*(q_s[1]*q_s[3] - q_s[0]*q_s[2])}
            //     }
            // };

            // // vector<vector<double> > y = la.add(b_list[0], v);
            // // vector<vector<double> > y_s = b_slist[0];
            // // z = la.add(y, la.mult(-1, y_s));
            
            z = la.add(la.dot(H, x), v);
            
            // update Kalman gain
            // vector<vector<double> > p = la.dot(la.T(H), la.dot(la.inv(R), H));
            // P = la.inv(la.add(la.inv(covM), p));
            // vector<vector<double> > sp1, sp2, p1, p2;
            // sp1 = la.dot(covM, la.T(H));
            // sp2 = la.dot(H, la.dot(covM, la.T(H)));
            // p1 = la.inv(la.add(sp2, R));
            // p2 = la.dot(sp1, la.dot(p1, la.dot(H, covM)));
            // P = la.add(covM, la.mult(-1, p2));
            // K = la.dot(P, la.dot(la.T(H), la.inv(R)));
            
            K = la.dot(covM, la.dot(la.T(H), la.inv(la.add(la.dot(H, la.dot(covM, la.T(H))), R))));
            P = la.dot(la.add(la.I(7), la.mult(-1, la.dot(K, H))), covM);
            
            // update x_hat
            vector<vector<double> > error = la.add(z, la.mult(-1, la.dot(H, x_bar)));
            x_hat = la.add(x_bar, la.dot(K, error));
            // x_hat = la.dot(K, z);

            // update estimation by system
            for(int i = 0; i < 7; i++){
                if(i <= 3) q_s[i] += x_hat[i][0];
                else omega_s[i - 3] += x_hat[i][0];
                x_hat[i][0] = 0;
            }
            q_s = d.normalize_q(la, q_s);
        }
        t += dt;
        idx++;
    }
    writing_file1.close();
    writing_file2.close();
    writing_file3.close();
    writing_file4.close();
    writing_file5.close();
    writing_file6.close();
    return 0;
}