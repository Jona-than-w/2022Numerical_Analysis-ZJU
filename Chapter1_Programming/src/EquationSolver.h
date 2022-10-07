#ifndef EQUATIONSOLVER_H
#define EQUATIONSOLVER_H

#include <iostream>
#include <cmath>
#include <limits>

//为了在diff方法利用差商法求导数
const double _epsilon = 20 * std::numeric_limits<double>::epsilon();

class Function{
public:
    virtual double operator () (const double &x) const = 0;
    double diff (const double &x) const{
        return ((*this)(x+_epsilon)-(*this)(x-_epsilon)) / (2*_epsilon);
    }
};

class EquationSolver{
public:
    virtual double solve() = 0;
};

class BisectionSolver : public EquationSolver{
private:
    double a, b, delta, eps;
    Function & f;
    int M;

public:
    BisectionSolver(Function & f, double a, double b, double delta, double eps, int M):
        f(f), a(a), b(b), delta(delta), eps(eps), M(M) {}
    double solve(){
        double u = f(a), v = f(b);
        double h = 0, c = 0, w = 0;
        for(int k = 1; k <= M; k ++){
            h = b-a;
            c = a + h/2;
            w = f(c);
            if(fabs(h) < delta || fabs(w) < eps)  
                break;
            else if(u * w < 0){
                b = c;
                v = w;
            }
            else{
                a = c;
                u = w;
            }
        }
        return c;
    }
};

class NewtonSolver : public EquationSolver{
private:
    double x0, eps;
    Function & f;
    int M;

public:
    NewtonSolver(Function & f, double x0, double eps, int M):
        f(f), x0(x0), eps(eps), M(M) {}
    double solve(){
        double x = x0, u, v;
        for(int k = 0; k <= M; k ++){
            u = f(x);
            if( fabs(u) < eps ) 
                break;
            v = u / f.diff(x);
            x = x - v;
        }
        return x;
    }
};

class SecantSolver : public EquationSolver{
private:
    double x0, x1, delta, eps;
    Function & f;
    int M;

public:
    SecantSolver(Function & f, double x0, double x1, double delta, double eps, int M):
        f(f), x0(x0), x1(x1), delta(delta), eps(eps), M(M) {}
    double solve(){
        int k;
        double u = f(x1), v = f(x0), s;
        for(k = 2; k <= M; k++){
            if( fabs(u) > fabs(v) ){
                std::swap(x0, x1);
                std::swap(u, v);
            }
            s = (x1 - x0) / (u - v);
            x0 = x1;
            v = u;
            x1 = x1 - u * s;
            u = f(x1);
            if( fabs(x1 - x0) < delta || fabs(u) < eps )
                break;
        }
        return x1;
    }
};

#endif
