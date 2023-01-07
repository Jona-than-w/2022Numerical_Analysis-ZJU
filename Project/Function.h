#ifndef _FUNCTION_H_
#define _FUNCTION_H_

#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

//ProblemA
/*class Function
{
    public:
    Function(){}
    vector<double> operator()(vector<double> _x){
        vector<double> _y;
        int N = _x.size();
        for (int i = 0; i < N; i++)
        {
            _y.push_back(1/(1+pow(_x[i],2)*25));
        }
        return _y;
    }
    double operator()(double _x)
    {
        return 1/(1+pow(_x,2)*25);
    }
    double df(double _x)
    {
        return -50*_x/(pow(1+pow(_x,2)*25,2));
    }
    double ddf(double _x)
    {
        return (3750*pow(_x,2)-50)/pow(1+pow(_x,2)*25,3);
    }
};*/


//ProblemBCD
class Function
{
    public:
    Function(){}
    vector<double> operator()(vector<double> _x){
        vector<double> _y;
        int N = _x.size();
        for (int i = 0; i < N; i++)
        {
            _y.push_back(1/(1+pow(_x[i],2)));
        }
        return _y;
    }
    double operator()(double _x)
    {
        return 1/(1+pow(_x,2));
    }
    double df(double _x)
    {
        return -2*_x/(pow(1+pow(_x,2),2));
    }
    double ddf(double _x)
    {
        return -2*(1-pow(_x,2))/pow(1+pow(_x,2),3);
    }
};

//ProblemE

/*class Function
{
public:
    Function() {}
    vector<double> operator()(vector<double> _x)
    {
        vector<double> _y;
        int N = _x.size();
        for (int i = 0; i < N; i++)
        {
   
            _y.push_back(3*cos(_x[i]));
        }
        
        return _y;
    }
    double operator()(double _x)
    {
        return (3*cos(_x));
    }
    double df(double _x)
    {
        return -3*sin(_x);
    }
    double ddf(double _x)
    {
        return -3*(cos(_x));
    }
};*/

#endif