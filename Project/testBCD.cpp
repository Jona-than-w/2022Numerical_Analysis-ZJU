#include"BSpline.h"
#include"pp.h"
#include"Function.h"
#include"Polynomial.h"

int main(){
    Function f1;

    int N;
    cin >> N;
    double a, b;
    cin >> a;
    cin >> b;

    vector<double> x, y;
    double h = (b - a) / (N - 1);
    for (int i = 0; i < N; i++)
    {
        x.push_back(a + h * i);
    }
    y = f1(x);

    BSpline s1(x, y,3);
    s1.complete(f1.df(x[0]), f1.df(x[N - 1]));
    // s2.specified_sec_diff(f1.ddf(x[0]),f1.ddf(x[N-1]),3);
    s1.Draw_Spline_with_Matlab_Code(N);
    s1.Draw_Error_with_Matlab_Code(N);
    cout << s1.error(-3.5) << " ";
    cout << s1.error(-3) << " ";
    cout << s1.error(-0.5) << " ";
    cout << s1.error(0) << " ";
    cout << s1.error(0.5) << " ";
    cout << s1.error(3) << " ";
    cout << s1.error(3.5) << " ";
    cout << endl;

    return 0;
}