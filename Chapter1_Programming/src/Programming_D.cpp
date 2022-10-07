#include "EquationSolver.h"
#include <limits>

const double pi = acos(-1);
const double eps = std::numeric_limits<double>::epsilon();

//----------------------------------------------------

class Func1 : public Function{
public:
    double operator () (const double &x) const{
        return sin(x/2) - 1.0;
    }
} func1;

double test1(){
    SecantSolver solver(func1, 0, pi/2, eps, eps, 100);
    double ans = solver.solve();
    std::cout << "The root is " << ans << " and the value of the function is " << func1(ans) << std::endl;
    return 0;
}

//----------------------------------------------------

class Func2 : public Function{
public:
    double operator () (const double &x) const{
        return exp(x) - tan(x);
    }
} func2;

double test2(){
    SecantSolver solver(func2, 1.0, 1.4, eps, eps, 100);
    double ans = solver.solve();
    std::cout << "The root is " << ans << " and the value of the function is " << func2(ans) << std::endl;
    return 0;
}

//----------------------------------------------------

class Func3 : public Function{
public:
    double operator () (const double &x) const{
        return x*x*x - 12*x*x + 3*x + 1;
    }
} func3;

double test3(){
    SecantSolver solver(func3, 0, -0.5, eps, eps, 100);
    double ans = solver.solve();
    std::cout << "The root is " << ans << " and the value of the function is " << func3(ans) << std::endl;
    return 0;
}

int main(){
    std::cout << "sin(x/2) - 1.0" << std::endl;
    test1();
    std::cout << "exp(x) - tan(x)" << std::endl;
    test2();
    std::cout << "x*x*x - 12*x*x + 3*x + 1" << std::endl;
    test3();
    return 0;
}