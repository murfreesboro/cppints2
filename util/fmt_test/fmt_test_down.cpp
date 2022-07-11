//
// This program is used to study how to calculate (SS|SS)^{m}
// for the recursive relation. Basically, it focus on the 
// fmt function:
// f_{m}(t) = \int^{1}_{0} u^{2m} e^{-tu^{2}} du 
// it needs the boost library to calcualte fmt
// fenglai liu
// Oct. 2013
//
// for reference, see Harris paper and the reference cited inside:
// Harris, Frank E
// Evaluation of GTO molecular integrals
// International Journal of Quantum Chemistry, 1983, Vol. 23, Page 1469--1478
//
// note:
// we found that the approximation of (2m-1)!!/R expression is 
// not numerical stable. Combined with down recursive relation,
// it introduces greater error compared with other methods. 
// Therefore we finally abadon it.
//
//
//

// common C head files used in the program
#include<cstdlib>
#include<cstdio>

// external math functions 
#include<cmath>

// common C++ head files may be used in the program
#include<iostream>
#include<fstream>
#include<string> 
#include<algorithm>

#include <boost/math/special_functions/binomial.hpp>  // special functions in math
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/timer.hpp>
using namespace boost::math;
using namespace boost;
using namespace std;
#define PI              3.1415926535897932384626E0  

long double fm(const long double& a, const int& m) {

	// if a is very small, the e^{-ax^2} term turns to zero
	if (fabs(a)<1.0E-20) return (1.0E0/(2*m+1));

	// this is the way we use incomplete gamma in boost
	long double f12 = 1.0E0/2.0E0;
	long double p   = 1.0E0/(2.0E0*pow(a,m+f12));
	long double upperLimit = a;
	return tgamma_lower(m+f12,upperLimit)*p;
}

int main() 
{
	///////////////////////////////////////////////////
	// global settings
	// step length and number of claculation for 
	// time performance
	// global error is the error range
	// if the difference is less than the error, we 
	// think the error could be omitted
	///////////////////////////////////////////////////
	long double steplength  = 0.00000001;
	long double globalError = 1.0E-16;
	long double T_max_limit = 40.0E0;
	int top_M_limmit        = 40;

	///////////////////////////////////////////////////
	// testing the fmt from down recursive relation
	///////////////////////////////////////////////////
	cout << endl << endl;
	cout << "===========================================================" << endl;
	cout << "testing the down recursive relation for m = 1 to " << top_M_limmit << endl;
	cout << "this is for long double type variable" << endl;
	cout << "error range is " << globalError << endl;
	cout << "step length for T is " << steplength << endl;
	cout << "T is ranging from 0 " << " to " <<T_max_limit << endl;
	cout << "===========================================================" << endl;
	size_t nSteps = T_max_limit/steplength;
	for (int n=1; n<=top_M_limmit; n++) {
		for (size_t j = 0; j<nSteps; j++) {
			long double T = 0.0E0 + steplength*j;
			long double fmt= fm(T,n);
			long double et = exp(-T);
			long double e0 = fmt;
			long double t2 = 2.0*T;
			for(int i=n-1; i>=0; i--) {
				long double s = fm(T,i);
				long double e = (1.0/(2*i+1.0))*(t2*e0+et);
				long double d = fabs(s-e);
				if (d>globalError) {
					printf("for m=%d T is %-18.16Lf difference is %-18.16Lf\n", n, T, d);
				}
				e0 = e;
			}
		}
	}

	return 0;
}
