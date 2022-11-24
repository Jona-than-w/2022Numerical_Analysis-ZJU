#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

class FPN{
private:
	int beta;
	int p;
	int L;
	int U;

public:
	FPN(int _beta = 2, int _p = 3, int _L = -1, int _U = 1){
		beta = _beta;
		p = _p;
		L = _L;
		U = _U;
	}
	
	void FPN_out(){
		double UFL = 1, OFL = 1;
		OFL = pow(beta, U) * (beta - pow(beta, 1-p));
		UFL = pow(beta, L);
		cout << "UFL = " << UFL << ", OFL = " << OFL << endl;
	}

	void get_normal(){
		vector<double> all_num;
		for (int j = L; j <= U; j++){
			double mantissa = 0;
			for (int k = 0; k < beta; k++)
				for (int l = 0; l < beta; l++)
				{
					mantissa = (2 * (2 + k)) + l;
					mantissa *= pow(beta, -2);	
					all_num.push_back(mantissa * pow(2, j));
					all_num.push_back(- mantissa * pow(2, j));
				}
		}
		all_num.push_back(0);
		for (int i = 0; i < all_num.size(); i++)
			cout << 0 << " " << all_num[i] << ";";
		cout << endl;
        //for (int i = 0; i < all_num.size(); i++)
		//	cout << 0 << " " << all_num[i] << ";";
		//cout << endl;
        //For the input type in matlab
	}

	void get_subnormal(){
		vector<double> all_num;
		for (int i = 0; i < beta; i++)
			for (int k = 0; k < beta; k++)
			{
				double mantissa = pow(2, L);
				mantissa *= (1.0 * i / 2 + 1.0 * k / 4);
				if (i + k != 0)
				{
					all_num.push_back(mantissa);
					all_num.push_back(-mantissa);
				}
			}
		for (int i = 0; i < all_num.size(); i++)
			cout << 0 << " " << all_num[i] << ";";
		cout << endl;
        //for (int i = 0; i < all_num.size(); i++)
		//	cout << 0 << " " << all_num[i] << ";";
		//cout << endl;
        //For the input type in matlab
	}
};

int main(){
	FPN b;

	b.FPN_out();
	b.get_normal();
	b.get_subnormal();
	return 0;
}

