#include "EquationSolver.h"
#include <limits>


const double pi = acos(-1);
const double eps = std::numeric_limits<double>::epsilon();

class Func1 : public Function{
public:
    double operator () (const double &x) const{
        return 1/x-tan(x);
    }
} func1;

double test1(){
    BisectionSolver solver(func1, 0, pi/2, eps, eps, 100);
    double ans = solver.solve();
    std::cout << ans << "  " << func1(ans) << std::endl;
    return 0;
}

//------------------------------------------------------------

class Func2 : public Function{
public:
    double operator () (const double &x) const{
        return 1/x-pow(2,x);
    }
} func2;

double test2(){
    BisectionSolver solver(func2, 0, 1, eps, eps, 100);
    double ans = solver.solve();
    std::cout << ans << "  " << func2(ans) << std::endl;
    return 0;
}

//-------------------------------------------------------------

class Func3 : public Function{
public:
    double operator () (const double &x) const{
        return pow(2, -x) + exp(x) + 2*cos(x) - 6;
    }
} func3;

double test3(){
    BisectionSolver solver(func3, 1, 3, eps, eps, 100);
    double ans = solver.solve();
    std::cout << ans << "  " << func3(ans) << std::endl;
    return 0;
}

//--------------------------------------------------------------

class Func4 : public Function{
public:
    double operator () (const double &x) const{
        return (x*x*x + 4*x*x + 3*x + 5) / (2*x*x - 9*x*x + 18*x - 2);
    }
} func4;

double test4(){
    BisectionSolver solver(func4, 0, 4, eps, eps, 100);
    double ans = solver.solve();
    std::cout << ans << "  " << func4(ans) << std::endl;
    return 0;
}

int main(){
    std::cout << "1.0/x-tan(x)" <<std::endl;
    test1();
    std::cout << "1.0/x-pow(2.0,x)" << std::endl;
    test2();
    std::cout << "pow(2.0, -x) + exp(x) + 2.0*cos(x) - 6.0" << std::endl;
    test3();
    std::cout << "(x*x*x + 4*x*x + 3*x + 5) / (2*x*x - 9*x*x + 18*x - 2)" << std::endl;
    test4();
    return 0;
}

