#include"BSpline.h"
#include"pp.h"
#include"Function.h"
#include"Polynomial.h"

/*class FuncA: Function{
public:
 vector<double> operator()(vector<double> _x){
        vector<double> _y;
        int N = _x.size();
        for (int i = 0; i < N; i++)
        {
            _y.push_back(1/(1+pow(_x[i],2)*25));
        }
        return _y;
    }
    double operator()(double _x){
        return 1/(1+pow(_x,2)*25);
    }
    double df(double _x){
        return -50*_x/(pow(1+pow(_x,2)*25,2));
    }
    double ddf(double _x){
        return (3750*pow(_x,2)-50)/pow(1+pow(_x,2)*25,3);
    }
} f1;*/

using namespace Eigen;
using namespace std;

int main(){
    Function f1;
    double a = -1, b = 1;
    vector<double> x, y;
    int N = 0;
    cin >> N;

        double h = (b - a) / (N - 1);
        for (int i = 0; i < N; i++){
            x.push_back(-1.0 + h * i);
        }
        y = f1(x);
        ppform s1(x, y);
        s1.complete_cubic(f1.df(x[0]), f1.df(x[N - 1]));
        //s1.specified_sec_diff(f1.ddf(x[0]),f1.ddf(x[N - 1]));
        s1.Draw_Spline_with_Matlab_Code(N);
        s1.Draw_Error_with_Matlab_Code(N);
        BSpline s2(x, y,3);
        s2.complete(f1.df(x[0]), f1.df(x[N - 1]));
        // s2.specified_sec_diff(f1.ddf(x[0]),f1.ddf(x[N-1]));
        s2.Draw_Spline_with_Matlab_Code(N);
        s2.Draw_Error_with_Matlab_Code(N);
        return 0;
}
