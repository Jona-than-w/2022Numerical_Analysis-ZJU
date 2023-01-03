#ifndef __POLYNOMIAL_H__
#define __POLYNOMIAL_H__

#include<iostream>
using namespace std;  
#include<vector>
#include<cmath>

class Polynomial{
    public:
        double LeftPoint;
        double RightPoint;
        int degree;
        vector<double> coef;
        double value(double x);
};

double Polynomial::value(double x){
    double result = coef[0];
    for (int i = 1; i <= degree; i++){
        result += coef[i]*pow(x-LeftPoint,i);
    }
    return result;
}

#endif 
