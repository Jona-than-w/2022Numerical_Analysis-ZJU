#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <iostream>
#include <cmath>
#include <vector>
using namespace std;
class Function{
public:
    virtual double operator()(double _x) {return 0;}
    virtual double diff(double _x) {return 0;}
};

Function func;

class NewtonPolynomial{
private:
    vector<double> x;
    vector<double> fvalue;
    vector<vector<double>> coef;
    int n;

public:
    NewtonPolynomial(vector<double> x, vector<double> fvalue, int n) : x(x), fvalue(fvalue), n(n){
        if (x.size() == fvalue.size()) ;
    }
    void solve(){
        vector<double> vec;
        vec.push_back(fvalue[0]);
        coef.push_back(vec);
        for (int i = 1; i < n; i++)
        {
            vec.clear();
            vec.push_back(fvalue[i]);
            double tmp;
            for (int j = 1; j <= i; j++){
                tmp = (vec[j - 1] - coef[i - 1][j - 1]) / (x[i] - x[i - j]);
                vec.push_back(tmp);
            }
            coef.push_back(vec);
        }
    }

    std::vector<std::vector<double>> outputcoef(){
        this->solve();
        return coef;
    }

    void outputlist(){
        this->solve();
        for (int i = 0; i < n; i++)
        {
            cout << coef[i][i] << endl;
        }
    }

    void outputtable()
    {
        this->solve();
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                cout << coef[i][j] << " ";
            }
            cout << "" << endl;
        }
    }

    void outputpoly(){
        this->solve();
        printf("%.5f",coef[0][0]);
        for (int i = 1; i < n; i++){
            if (coef[i][i] > 0)
                printf("+%.5f", coef[i][i]);
            else
                printf("%.5f", coef[i][i]);
            // cout << coef[i][i];
            for (int j = 0; j < i; j++)
            {
                if (x[j] > 0)
                    // printf("(x-%f)", x[j]);
                    cout << "(x" << -x[j] << ")";
                else if (x[j] < 0)
                    // printf("(x+%f)", -x[j]);
                    cout << "(x+" << -x[j] << ")";
                else
                    cout << "x";
            }
        }
    }
};

int factorial(int _k)
{
    int tmp = 1;
    if (_k == 0)
        return tmp;
    else
    {
        for (int i = 1; i <= _k; i++)
            tmp = tmp * i;
    }
    return tmp;
}


class HermitePolynomial
{
private:
    vector<double> x;
    vector<double> fvalue;
    vector<vector<double>> coef;
    int n;

public:
    HermitePolynomial(vector<double> _x, vector<double> _fvalue, int _n) : x(_x), fvalue(_fvalue), n(_n) {}

    void solve()
    {
        vector<double> vec;
        int cnt = 0;
        vec.push_back(fvalue[0]);
        coef.push_back(vec);
        for (int i = 1; i < n; i++)
        {
            if (x[i] == x[i - 1])
                cnt++;
            else
                cnt = 0;
            vec.clear();
            double tmp;
            for (int k = 0; k <= cnt; k++)
            {
                vec.push_back(fvalue[i - cnt + k] / factorial(k));
            }
            for (int j = cnt + 1; j <= i; j++)
            {
                double y1 = coef[i - 1][j - 1];
                double y2 = vec[j - 1];
                tmp = (y2 - y1) / (x[i] - x[i - j]);
                vec.push_back(tmp);
            }
            coef.push_back(vec);
        }    
    }
    
    std::vector<vector<double>> outputcoef(){
        this->solve();
        return coef;
    }

    void outputlist(){
        this->solve();
        for (int i = 0; i < n; i++)
        {
            cout << coef[i][i] << endl;
        }
    }

    void outputtable(){
        this->solve();
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                cout << coef[i][j] << " ";
            }
            cout << "" << endl;
        }
    }

    void outputpoly(){
        this->solve();
        cout << coef[0][0];
        for (int i = 1; i < n; i++)
        {
            if (coef[i][i] >= 0)
                printf("+%.5f", coef[i][i]);
            // cout << "+" << coef[i][i];
            else
                printf("%.5f", coef[i][i]);
            // cout << coef[i][i];
            for (int j = 0; j < i; j++)
            {
                if (x[j] > 0)
                    // printf("*(x-%d)", x[j]);
                    cout << "(x" << -x[j] << ")";
                else if (x[j] < 0)
                    // printf("*(x+%d)", -x[j]);
                    cout << "(x+" << -x[j] << ")";
                else
                    cout << "*x";
            }
        }
    }
};

#endif