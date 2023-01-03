#ifndef _FUNCTION_H_
#define _FUNCTION_H_

#include <limits>
#include <iostream>
#include <cstring>
#include <cmath>

using namespace std;
const double _epsilon = 10 * std::numeric_limits<double>::epsilon();

class Function{
public:
    virtual vector<double> operator () (vector<double> _x){return {0};}
    virtual double operator () (const double &x){return 0;}

    virtual double diff(double x){
        return ((*this)(x + _epsilon) - (*this)(x - _epsilon)) / (2 * _epsilon);
    }
    
    virtual double second_diff(double x){
        return ((*this)(x + _epsilon) - 2*(*this)(x) + (*this)(x - _epsilon)) / (2 * _epsilon * _epsilon);
    }
    
};

#endif