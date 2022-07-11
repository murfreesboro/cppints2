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
// july 2022:
// double check the implementation and everything good, revise the fmt function
// and use long double instead.
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
#include<vector>
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
	int M_limmit = 4;
	int M_lower_limit = 1;
	long double T_min_limit = 1.39E0;
   int nTerms = 17;

	// testing the accurary of long double
	printf("sizeof long double: %d\n", sizeof(long double));

	///////////////////////////////////////////////////
	// testing the power series
	// equation of 9 in Harris paper
	// We try to expand the power series 
	///////////////////////////////////////////////////
	size_t numSteps = T_min_limit/steplength;
	cout << endl << endl;
	cout << "===========================================================" << endl;
	cout << "testing the power series expansion with " << nTerms+1 << " terms" << endl;
	cout << "this is for long double type variable" << endl;
	cout << "error range is " << globalError << endl;
	cout << "step length for T is " << steplength << endl;
	cout << "T is ranging from 0 to " << T_min_limit << endl;
	cout << "m is from "  << M_lower_limit << " to " << M_limmit << endl;
	cout << "===========================================================" << endl;
	for (int m = M_lower_limit; m<=M_limmit; m++) {
		for (size_t j = 0; j<=numSteps; j++) {
			long double T = 0.0E0 + steplength*j;
			long double et = exp(-T);
			long double x  = 1.0;
			long double t2 = 2.0*T;
			for(int i=nTerms; i>=1; i--) {
				x = 1.0+t2/(2.0*(m+i)+1.0)*x;
			}
			x  = (1.0/(2*m+1))*x;
			x *= et;
			long double result = fm(T,m);
			long double d = fabs(result-x);
			if (d>globalError) {
				printf("fmt is %-18.16Lf, for m=%d T=%-18.16Lf difference is %-18.16Lf\n", result, m, T, d);
			}
		}
	}

	return 0;
}
