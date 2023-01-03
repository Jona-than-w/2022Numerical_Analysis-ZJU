#ifndef _SPLINE_H_
#define _SPLINE_H_

#include <iostream>
#include <fstream>
#include <math.h>
#include "Eigen/Dense"
#include <vector>

using namespace std;
using namespace Eigen;

class ppform{
    private:
        int N;
        vector<double> x, y, h; 
        vector<double> m, M, K;//s'(f,x), s''(f,x), f[x_i,x_{i+1}]
        vector<double> lambda, miu;
        MatrixXd A;

    public:
        ppform() {}

        ppform(vector<double> _x, vector<double> _y){
            if(_x.size() != _y.size()) std::cerr << "The input x doesn't match y in size!!" << std::endl;
            else{
                    //initializing
                x = _x;
                y = _y;
                N = x.size();
                h.push_back(x[1] - x[0]);
                //initialize the number characteristic
                for (int i = 1; i <= N - 2; i++){
                    h.push_back(x[i + 1] - x[i]);
                    lambda.push_back((x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]));
                    miu.push_back((x[i + 1] - x[i]) / (x[i + 1] - x[i - 1]));
                }

                A = MatrixXd::Zero(N - 2, N - 2);
                for (int i = 0; i <= N - 3; i++){
                    A(i, i) = 2;
                    A(i, i + 1) = miu[i];
                    A(i + 1, i) = lambda[i + 1];
                }

                for (int i = 0; i < N - 1; i++){
                    K.push_back((y[i + 1] - y[i]) / h[i]);
                }
            }
        }

        void complete_cubic(double m_1, double m_N){

            // set vector b from 2 to n-1
            VectorXd b(N - 2);

            for (int i = 0; i <= N - 3; i++){
                b(i) = 3.0 * (miu[i] * (y[i + 2] - y[i + 1]) / h[i + 1] + lambda[i] * (y[i + 1] - y[i]) / h[i]);
            }

            b(0) = b(0) - m_1 * lambda[0];
            b(N - 3) = b(N - 3) - m_N * miu[N - 3];

            // solve vector X
            VectorXd X = A.lu().solve(b);
            m.push_back(m_1);
            for (int i = 0; i < N - 2; i++)
                m.push_back(X(i));
            m.push_back(m_N);
        }

        void specified_sec_diff(double M_1, double M_N){
            VectorXd b(N - 2);
            for (int i = 0; i <= N - 3; i++)
                b(i) = 6.0 * ((y[i + 2] - y[i + 1]) / h[i + 1] - (y[i + 1] - y[i]) / h[i]) / (h[i + 1] + h[i]);
            b(0) = b(0) - miu[0] * M_1;
            b(N - 3) = b(N - 3) - lambda[N - 3] * M_N;
            VectorXd X = A.lu().solve(b);
            M.push_back(M_1);
            for (int i = 0; i < N - 2; i++)
                M.push_back(X(i));
            M.push_back(M_N);
            for (int i = 0; i < N - 1; i++)
            {
                m.push_back(K[i] - 1.0 / 6 * (M[i + 1] + 2 * M[i]) * h[i]);
            }
            m.push_back(K[N - 2] + 1.0 / 6 * (M[N - 2] + 2 * M[N - 1]) * h[N - 2]);
        }

        void natural_cubic(){
            specified_sec_diff(0.0, 0.0);
        }

        double cubic(int _i, double _x) {
            double a = (K[_i - 1] - m[_i - 1]) / h[_i - 1];
            double b = (m[_i - 1] + m[_i] - 2 * K[_i - 1]) / (h[_i - 1] * h[_i - 1]);
            return y[_i - 1] + m[_i - 1] * (_x - x[_i - 1]) + a * pow(_x - x[_i - 1], 2) + b * pow(_x - x[_i - 1], 2) * (_x - x[_i]);
        }

        void Draw_with_Matlab_Code()
        {
            ofstream outfile;
            outfile.open("Code_ppform.m");
            outfile << "%% Draw with ppform " << endl;
            outfile << "x= [";
            for (int i = 1; i <= N - 1; i++)
            {
                double delta = h[i - 1] / 500;
                for (int j = 0; j < 500; j++)
                    outfile << x[i - 1] + j * delta << ",";
            }
            outfile << x[N - 1] << "];" << endl;

            outfile << " y" << N << " = [";
            for (int i = 1; i <= N - 1; i++)
            {
                double delta = h[i - 1] / 500;
                for (int j = 0; j < 500; j++)
                    outfile << cubic(i, x[i - 1] + j * delta) << ",";
            }
            outfile << y[N - 1] << "];" << endl;
            outfile << "plot(x,y" << N << ");" << endl;
            outfile.close();
        }
};

#endif