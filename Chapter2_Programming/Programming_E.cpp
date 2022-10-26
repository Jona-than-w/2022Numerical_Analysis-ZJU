#include <cmath>
#include "interpolation.h"

int n = 7;
double xvalues[] = {0, 6, 10, 13, 17, 20, 28};
double fvalues[] = {6.67, 17.3, 42.7, 37.3, 30.1, 29.3, 28.7};
double fvalues1[] = {6.67, 16.1, 18.9, 15.0, 10.6, 9.44, 8.89};
vector<double> x(xvalues, xvalues + n);
vector<double> fvalue(fvalues, fvalues + n);
vector<double> fvalue1(fvalues1, fvalues1 + n);


int main(){
    HermitePolynomial Pfunc(x, fvalue, n);
    Pfunc.outputpoly();
    cout << endl;
    HermitePolynomial Pfunc1(x, fvalue1, n);
    Pfunc1.outputpoly();
    return 0;
}