#ifndef _PP_H_
#define _PP_H_

#include <iostream>
#include <fstream>
#include <math.h>
#include"Function.h"
#include "Eigen/Dense"
#include <vector>
#define CONSTN 500

using namespace std;
using namespace Eigen;
Function f;

class ppform
{
private:
    int N;
    vector<double> x, y, h, midpoints; // x-value, y-value, the interval
    vector<double> m, M, K;
    vector<double> lambda, mu;
    MatrixXd A;

public:
    ppform() {}
    ppform(vector<double> _x, vector<double> _y)
    {
        x = _x;
        y = _y;
        if (x.size() != y.size())
            cout << "Error!" << endl;
        N = x.size();
        h.push_back(x[1] - x[0]);
        midpoints.push_back(x[0] + h[0] / 2);
        for (int i = 1; i <= N - 2; i++)
        {
            h.push_back(x[i + 1] - x[i]);
            midpoints.push_back(x[i] + h[i] / 2);
            lambda.push_back((x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]));
            mu.push_back((x[i + 1] - x[i]) / (x[i + 1] - x[i - 1]));
        }
        // set the matrix A
        A = MatrixXd::Zero(N - 2, N - 2);
        for (int i = 0; i <= N - 3; i++)
        {
            A(i, i) = 2;
        }
        for (int i = 0; i < N - 3; i++)
        {
            A(i, i + 1) = mu[i];
            A(i + 1, i) = lambda[i + 1];
        }
        for (int i = 0; i < N - 1; i++)
        {
            K.push_back((y[i + 1] - y[i]) / h[i]);
        }
    }
    void complete_cubic(double m_1, double m_N)
    {

        // set vector b
        VectorXd b(N - 2);

        for (int i = 0; i <= N - 3; i++)
        {
            b(i) = 3.0 * (mu[i] * (y[i + 2] - y[i + 1]) / h[i + 1] + lambda[i] * (y[i + 1] - y[i]) / h[i]);
        }
        b(0) = b(0) - m_1 * lambda[0];
        b(N - 3) = b(N - 3) - m_N * mu[N - 3];
        // solve vector X
        VectorXd X = A.lu().solve(b);
        m.push_back(m_1);
        for (int i = 0; i < N - 2; i++)
            m.push_back(X(i));
        m.push_back(m_N);
    }

    void specified_sec_diff(double M_1, double M_N)
    {
        VectorXd b(N - 2);
        for (int i = 0; i <= N - 3; i++)
            b(i) = 6.0 * ((y[i + 2] - y[i + 1]) / h[i + 1] - (y[i + 1] - y[i]) / h[i]) / (h[i + 1] + h[i]);
        b(0) = b(0) - mu[0] * M_1;
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

    void natural_cubic()
    {
        specified_sec_diff(0.0, 0.0);
    }
    double cubic(int _i, double _x) // cubic(i,x)代表第i段，一共N-1段
    {
        double a = (K[_i - 1] - m[_i - 1]) / h[_i - 1];
        double b = (m[_i - 1] + m[_i] - 2 * K[_i - 1]) / (h[_i - 1] * h[_i - 1]);
        return y[_i - 1] + m[_i - 1] * (_x - x[_i - 1]) + a * pow(_x - x[_i - 1], 2) + b * pow(_x - x[_i - 1], 2) * (_x - x[_i]);
    }

    void Draw_Spline_with_Matlab_Code(int i)
    {
        ofstream outfile;
        string filename = "Code_ppform.m";//have to change the name of matlab
        string tmp = to_string(i);
        filename.insert(11,tmp);
        outfile.open(filename,ios::app);
        outfile << "%% Draw with ppformm " << endl;
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

        outfile << " y0"
                << " = [";
        for (int i = 1; i <= N - 1; i++)
        {
            double delta = h[i - 1] / 500;
            for (int j = 0; j < 500; j++)
                outfile << f(x[i - 1] + j * delta) << ",";
        }
        outfile << y[N - 1] << "];" << endl;
        outfile << "plot(x,y" << N << ",x,y0);" << endl;
        outfile.close();
    }
    void Draw_Error_with_Matlab_Code(int i)
    {
        ofstream outfile;
        string filename = "Code_ppform_error.m";//have to change the name of matlab
        string tmp = to_string(i);
        filename.insert(17,tmp);
        outfile.open(filename,ios::app);
        outfile << "%% Draw with ppformm " << endl;
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
            {
                outfile << abs(cubic(i, x[i - 1] + j * delta) - f(x[i - 1] + j * delta)) << ",";
            }
        }
        outfile << 0 << "];" << endl;
        outfile << "midpoints_error = [";
        for (int i = 0; i < N - 1; i++)
            outfile << abs(cubic(i + 1, midpoints[i]) - f(midpoints[i])) << ",";
        outfile << "0];" << endl;
        outfile << "m = max(midpoints_error);" << endl;
        outfile << "plot(x,y" << N << ");" << endl;
        outfile.close();
    }
    
    void linear()
    {
        ofstream outfile;
        outfile.open("Code_Linear.m", ios::app);
        outfile << "%% Draw Linear" << endl;
        outfile << "x= [";
        for (int i = 1; i <= N; i++)
        {
            outfile << x[i + 2] << ",";
        }
        outfile << "];" << endl;

        outfile << " y" << N << " = [";
        for (int i = 1; i <= N; i++)
        {
            outfile << f(x[i + 2]) << ",";
        }
        outfile << "];" << endl;
        outfile << "plot(x,y" << N << ");" << endl;
        outfile.close();
    }
};

#endif