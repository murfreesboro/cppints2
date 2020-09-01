#include "general.h"
#include "basis.h"
#include "shell.h"
#include "integral.h"
#include "shellquartet.h"
using namespace integral;
using namespace basis;
using namespace shell;
using namespace shellquartet;

int main() {

	Basis test(0,2,0);
	cout << test.getName() << endl;
	int oper = ERI;
	Integral I1(test,test,test,test,oper,100);
	cout << I1.getName() << endl;

	Shell s(3);
	vector<Basis> b;
	s.getBasis(b);
	for(int i=0; i<(int)b.size(); i++)
		cout << b[i].getName() << endl;

	if (s.hasThisBasisSet(test)) {
		cout << "yes!!" << endl;
	}

	ShellQuartet sq1(2,2,2,2,ERI,0);
	cout << sq1.getName() << endl;
	ShellQuartet sq2(2,1,2,1,ERI,0);
	cout << sq2.getName() << endl;
	if (sq1<sq2) {
		cout << "less than!!!!" << endl;
	}else{
		cout << "larger than!!!" << endl;
	}

	Integral I2(sq1,6);
	cout << I2.getName() << endl;

}
