#include <cmath>
#include "interpolation.h"

int n = 10;
double xvalues[] = {0, 0, 3, 3, 5, 5, 8, 8, 13, 13};
double fvalues[] = {0, 75, 225, 77, 383, 80, 623, 74, 993, 72};
vector<double> x(xvalues, xvalues + n);
vector<double> fvalue(fvalues, fvalues + n);

int main(){
    HermitePolynomial Pfunc(x, fvalue, n);
    Pfunc.outputpoly();
    return 0;
}