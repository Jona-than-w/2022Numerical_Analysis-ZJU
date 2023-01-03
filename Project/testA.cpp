#include"BSpline.h"
#include"pp.h"
#include"Function.h"
#include"Polynomial.h"

class FuncA: Function{
public:
    vector<double> operator () (const vector<double> _x) const{
        vector<double> x = _x;
        vector<double> _y;
        int N = _x.size();
        for (int i = 0; i < N; i++){
            _y.push_back(1.0/(1.0+25*x[i]*x[i]));
        }
        return _y;
    }
    double operator () (const double &x) const{
        return 1.0/(1.0+25*x*x);
    }
    double diff (const double &x) const{
        return -50.0*x/pow(1.0+25*x*x,2);
    }

    double second_diff(double _x){
        return (3750*pow(_x,2)-50)/pow(1+pow(_x,2)*25,3);
    }
} f1;

int main(){
    double a = -1, b = 1;
    vector<double> x, y;

    int N[5] = {6,11,21,41,81};
    for(int j = 0; j < 5; j ++){
        double h = (b - a) / (N[j] - 1);
        for (int i = 0; i < N[j]; i++){
            x.push_back(-1.0 + h * i);
        }
        y = f1(x);
        ppform s1(x, y);
        s1.complete_cubic(f1.diff(x[0]), f1.diff(x[N[j] - 1]));
        s1.specified_sec_diff(f1.second_diff(x[0]),f1.second_diff(x[N[j] - 1]));
        s1.Draw_Spline_with_Matlab_Code();
        s1.Draw_Error_with_Matlab_Code();
    }
}