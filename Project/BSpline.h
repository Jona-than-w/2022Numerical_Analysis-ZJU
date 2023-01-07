#ifndef _BPLINE_H_
#define _BPLINE_H_

#include <iostream>
#include <fstream>
#include <math.h>
#include "Eigen/Dense"
#include <vector>
#include "Function.h"

using namespace std;
using namespace Eigen;

Function _f;

class BSpline
{
private:
    int N;
    int k;
    vector<double> x, y, h, midpoints; // x-value, y-value, the interval
    vector<double> coef;

public:
    BSpline() {}
    BSpline(vector<double> _x, vector<double> _y, int _k)
    {
        x = _x;
        y = _y;
        k = _k;
        if (x.size() != y.size())
            cout << "Error!" << endl;
        N = _x.size();
        h.push_back(x[1] - x[0]);
        midpoints.push_back(x[0] + h[0] / 2);
        for (int i = 1; i <= N - 2; i++)
        {
            h.push_back(x[i + 1] - x[i]);
            midpoints.push_back(x[i] + h[i] / 2);
        }
        for (int i = 1; i <= 3; i++)
        {
            x.insert(x.begin(), _x[0] - i);
            x.insert(x.end(), _x[N - 1] + i);
        }
    }
    double B(int _i, int _n, double _x) // B_i^n(x)定义和书上一致 B_1^0(x)
    {

        if (_n == 0)
        {

            if (_x > x[_i + 1] && _x <= x[_i + 2])
                return 1;
            else
                return 0;
        }
        else
        {
            return ((_x - x[_i + 1]) / (x[_i + _n + 1] - x[_i + 1]) * B(_i, _n - 1, _x) + (x[_i + _n + 2] - _x) / (x[_i + _n + 2] - x[_i + 2]) * B(_i + 1, _n - 1, _x));
        }
    }

    double dB(int _i, int _n, double _x)
    {
        if (_n == 0)
        {
            return 0;
        }
        else
        {
            return (1 / (x[_i + _n + 1] - x[_i + 1]) * B(_i, _n - 1, _x) + (_x - x[_i + 1]) / (x[_i + _n + 1] - x[_i + 1]) * dB(_i, _n - 1, _x) + (x[_i + _n + 2] - _x) / (x[_i + _n + 2] - x[_i + 2]) * dB(_i + 1, _n - 1, _x) - 1 / (x[_i + _n + 2] - x[_i + 2]) * B(_i + 1, _n - 1, _x));
        }
    }

    double d2B(int _i, int _n, double _x)
    {
        if (_n == 0)
        {
            return 0;
        }
        else
        {
            return (2.0 / (x[_i + _n + 1] - x[_i + 1]) * dB(_i, _n - 1, _x) + (_x - x[_i + 1]) / (x[_i + _n + 1] - x[_i + 1]) * d2B(_i, _n - 1, _x) + (x[_i + _n + 2] - _x) / (x[_i + _n + 2] - x[_i + 2]) * d2B(_i + 1, _n - 1, _x) - 2.0 / (x[_i + _n + 2] - x[_i + 2]) * dB(_i + 1, _n - 1, _x));
        }
    }

    void complete(double m_1, double m_N)
    {
        MatrixXd A = MatrixXd::Zero(N + 2, N + 2);
        for (int i = 0; i < N; i++)
        {
            A(i, i) = B(i - 1, k, x[i + 3]);
            A(i, i + 1) = B(i, k, x[i + 3]);
            A(i, i + 2) = B(i + 1, k, x[i + 3]);
        }
        A(N, 0) = dB(-1, k, x[3]);
        A(N, 1) = dB(0, k, x[3]);
        A(N, 2) = dB(1, k, x[3]);
        A(N + 1, N - 1) = dB(N - 2, k, x[N + 2]);
        A(N + 1, N) = dB(N - 1, k, x[N + 2]);
        A(N + 1, N + 1) = dB(N, k, x[N + 2]);

        VectorXd b(N + 2);
        for (int i = 0; i < N; i++)
            b(i) = y[i];
        b(N) = m_1;
        b(N + 1) = m_N;
        VectorXd a = A.lu().solve(b);
        for (int i = 0; i < N + 2; i++)
        {
            coef.push_back(a(i));
        }
    }

    void specified_sec_diff(double M_1, double M_N)
    {
        MatrixXd A = MatrixXd::Zero(N + 2, N + 2);
        for (int i = 0; i < N; i++)
        {
            A(i, i) = B(i - 1, k, x[i + 3]);
            A(i, i + 1) = B(i, k, x[i + 3]);
            A(i, i + 2) = B(i + 1, k, x[i + 3]);
        }
        A(N, 0) = d2B(-1, k, x[3]);
        A(N, 1) = d2B(0, k, x[3]);
        A(N, 2) = d2B(1, k, x[3]);
        A(N + 1, N - 1) = d2B(N - 2, k, x[N + 2]);
        A(N + 1, N) = d2B(N - 1, k, x[N + 2]);
        A(N + 1, N + 1) = d2B(N, k, x[N + 2]);

        VectorXd b(N);
        for (int i = 0; i < N; i++)
            b(i) = y[i];
        b(N) = M_1;
        b(N + 1) = M_N;
        VectorXd a = A.lu().solve(b);
        for (int i = 0; i < N + 2; i++)
            coef.push_back(a(i));
    }

    void natural_cubic()
    {
        specified_sec_diff(0.0, 0.0);
    }

    double compute(double _x)
    {
        double _y = 0;
        for (int i = -1; i <= N; i++)
        {
            _y = _y + coef[i + 1] * B(i, k, _x);
        }
        return _y;
    }

    double error(double _x)
    {
        return abs(compute(_x) - _f(_x));
    }

    void Draw_Spline_with_Matlab_Code(int i)
    {
        ofstream outfile;
        string filename = "Code_Bspline.m";//have to change the name of matlab
        string tmp = to_string(i);
        filename.insert(12,tmp);
        outfile.open(filename,ios::app);
        outfile << "%% Draw with BSpline" << endl;
        outfile << "x= [";
        for (int i = 1; i <= N - 1; i++)
        {
            double delta = h[i - 1] / 500;
            for (int j = 0; j < 500; j++)
                outfile << x[i + 2] + j * delta << ",";
        }
        outfile << x[N + 2] << "];" << endl;

        outfile << " y" << N << " = [";
        for (int i = 1; i <= N - 1; i++)
        {
            double delta = h[i - 1] / 500;
            for (int j = 0; j < 500; j++)
                outfile << compute(x[i + 2] + j * delta) << ",";
        }
        outfile << y[N - 1] << "];" << endl;
        outfile << " y0"
                << " = [";
        for (int i = 1; i <= N - 1; i++)
        {
            double delta = h[i - 1] / 500;
            for (int j = 0; j < 500; j++)
                outfile << _f(x[i + 2] + j * delta) << ",";
        }
        outfile << y[N - 1] << "];" << endl;
        outfile << "plot(x,y" << N << ",x,y0);" << endl;
        outfile.close();
    }

    void Draw_Error_with_Matlab_Code(int i)
    {   
        ofstream outfile;
        string filename = "Code_Bspline_error.m";//have to change the name of matlab
        string tmp = to_string(i);
        filename.insert(18,tmp);
        outfile.open(filename,ios::app);
        outfile << "%% Draw with BSpline Error" << endl;
        outfile << "x= [";
        for (int i = 1; i <= N - 1; i++)
        {
            double delta = h[i - 1] / 500;
            for (int j = 0; j < 500; j++)
                outfile << x[i + 2] + j * delta << ",";
        }
        outfile << x[N + 2] << "];" << endl;

        outfile << " y" << N << " = [";
        for (int i = 1; i <= N - 1; i++)
        {
            double delta = h[i - 1] / 500;
            for (int j = 0; j < 500; j++)
                outfile << abs(compute(x[i + 2] + j * delta) - _f(x[i + 2] + j * delta)) << ",";
        }
        outfile << 0 << "];" << endl;

        outfile << "midpoints_error = [";
        for (int i = 0; i < N - 1; i++)
            outfile << abs(compute(midpoints[i]) - _f(midpoints[i])) << ",";
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
            outfile << _f(x[i + 2]) << ",";
        }
        outfile << "];" << endl;
        outfile << "plot(x,y" << N << ");" << endl;
        outfile.close();
    }
};

#endif
