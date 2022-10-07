#include "EquationSolver.h"
#include <limits>

const double eps = std::numeric_limits<double>::epsilon();

class Func : public Function{
public:
    double operator () (const double &x) const{
        return x-tan(x);
    }
    double diff (const double &x) const{
        return 1 - 1.0/(cos(x)*cos(x));//用了数学意义上的导函数而非差商计算
    }
} f;

int main(){
    NewtonSolver solver1(f, 4.5, eps, 100);
    double ans1 = solver1.solve();
    NewtonSolver solver2(f, 7.7, eps, 100);
    double ans2 = solver2.solve();
    std::cout << "The root near 4.5 is " << ans1 << std::endl;
    std::cout << "The root near 7.7 is " << ans2 << std::endl;
    return 0;
}