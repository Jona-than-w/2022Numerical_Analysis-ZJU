#include <cmath>
#include "interpolation.h"
class Func : public Function{
    public:
    double operator()(double x)
    {
        return 1 / (1 + x * x);
    }
} func0;

vector<double> x, fvalue;
const int n = 8;

int main(){
    for (int i = 0; i <= n; i++){
        x.push_back(-5 + 10.0 * i / n);
        fvalue.push_back(func0(x[i]));
    }

    NewtonPolynomial Pfunc(x, fvalue, n + 1);
    //Pfunc.outputlist();
    //cout << endl;
    //Pfunc.outputtable();
    //cout << endl;
    Pfunc.outputpoly(); 
    std::cout << endl;
    std::cout << endl;
}
