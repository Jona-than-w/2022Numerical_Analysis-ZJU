#include"BSpline.h"
#include"pp.h"
#include"Function.h"
#include"Polynomial.h"


int main()
{
    Function f1;

    int N ;
    cin >> N;
    double pi= acos(-1);
    double a =-pi/2, b = 3*pi/2;
    vector<double> x, y,theta;
    double h = (b - a) / (N - 1);
    for (int i = 0; i < N; i++)
    {
        x.push_back(a + h * i);
    }
    y = f1(x);
    ppform s1(x, y);
    //s1.complete_cubic(f1.df(x[0]), f1.df(x[N - 1]));
    s1.natural_cubic();
    s1.Draw_Spline_with_Matlab_Code(N);

    return 0;
}