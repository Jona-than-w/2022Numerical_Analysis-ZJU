#include <cmath>
#include "interpolation.h"
double pi = acos(-1);
class Func : public Function{
    public:
    double operator()(double x)
    {
        return 1 / (1 + 25 * x * x);
    }
}func0;

vector<double> x, fvalue;
int n = 15;

int main(){
   for (int i = 1; i <= n; i++)
    {
        x.push_back(cos((2 * i - 1) / (2.0 * n) * pi));
        fvalue.push_back(func0(x[i - 1]));
    }
    HermitePolynomial Pfunc(x, fvalue, n);
    Pfunc.outputpoly();
    return 0;
}