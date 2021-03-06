//
// CPPINTS: A C++ Program to Generate Analytical Integrals Based on Gaussian
// Form Primitive Functions
//
// Copyright (C) 2015 The State University of New York at Buffalo
// This softare uses the MIT license as below:
//
//	Permission is hereby granted, free of charge, to any person obtaining 
//	a copy of this software and associated documentation files (the "Software"), 
//	to deal in the Software without restriction, including without limitation 
//	the rights to use, copy, modify, merge, publish, distribute, sublicense, 
//	and/or sell copies of the Software, and to permit persons to whom the Software 
//	is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all 
// copies or substantial portions of the Software.
//						    
//	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
//	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
//	PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
//	FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
//	ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//
#include "shell.h"
#include "inttype.h"
#include "boost/lexical_cast.hpp"
#include "shellquartet.h"
#include "derivinfor.h"
#include "integral.h"
using namespace inttype;
using namespace shell;
using boost::lexical_cast;
using namespace shellquartet;
using namespace derivinfor;
using namespace integral;

Integral::Integral(const ShellQuartet& sq, const int& index)
{

	// determine the basis set length
	const Shell& sbra1 = sq.getShell(BRA1); 
	const Shell& sbra2 = sq.getShell(BRA2); 
	const Shell& sket1 = sq.getShell(KET1); 
	const Shell& sket2 = sq.getShell(KET2); 
	int n1 = 0;
	int n2 = 0;
	int n3 = 0;
	int n4 = 0;
	n1 = sbra1.getBasisSetNumber();
	if (! sbra2.isnull()) n2 = sbra2.getBasisSetNumber();
	if (! sket1.isnull()) n3 = sket1.getBasisSetNumber();
	if (! sket2.isnull()) n4 = sket2.getBasisSetNumber();

	// we try to get the sub index of basis set in shell
	int a1 = 0;  // bra1
	int a2 = 0;  // bra2
	int a3 = 0;  // ket1
	int a4 = 0;  // ket2
	if (index > 0) {
		if (index<n1) {
			a1 = index;
		}else if (index>=n1 && index<n1*n2) {
			a1 = index%n1;
			a2 = (index - a1)/n1;
		}else if (index>=n1*n2 && index<n1*n2*n3) {
			int nt = n1*n2;
			int tmp = index%nt;
			a3 = (index-tmp)/nt;
			a1 = tmp%n1;
			a2 = (tmp - a1)/n1;
		}else if (index >= n1*n2*n3 && index < n1*n2*n3*n4){

			// get the a4
			int nt1 = n1*n2*n3;
			int t1  = index%nt1;
			a4 = (index-t1)/nt1;

			// get the a3, a2 and a1
			int nt2 = n1*n2;
			int t2  = t1%nt2;
			a3 = (t1-t2)/nt2;
			a1 = t2%n1;
			a2 = (t2 - a1)/n1;
		}else{
			crash(true,"Wrong index can not be parsed in integral constructor");
		}
	}
	//cout << a1 << " " << a2 << " " << a3 << " " << a4 << endl;

	// now let's get the basis set
	int l = -1;
	int m = -1;
	int n = -1;
	sbra1.getBasisSetFromIndex(a1,l,m,n);
	Basis tbra1(l,m,n);
	bra1 = tbra1;
	if (n2 > 0) { 
		sbra2.getBasisSetFromIndex(a2,l,m,n);
		Basis tbra2(l,m,n);
		bra2 = tbra2;
	}else{
		Basis tbra2(NULL_POS,NULL_POS,NULL_POS);
		bra2 = tbra2;
	}
	if (n3 > 0) { 
		sket1.getBasisSetFromIndex(a3,l,m,n);
		Basis tket1(l,m,n);
		ket1 = tket1;
	}else{
		Basis tket1(NULL_POS,NULL_POS,NULL_POS);
		ket1 = tket1;
	}
	if (n4 > 0) {
		sket2.getBasisSetFromIndex(a4,l,m,n);
		Basis tket2(l,m,n);
		ket2 = tket2;
	}else{
		Basis tket2(NULL_POS,NULL_POS,NULL_POS);
		ket2 = tket2;
	}
	//cout << bra1.getName() << endl;
	//cout << bra2.getName() << endl;
	//cout << ket1.getName() << endl;
	//cout << ket2.getName() << endl;

	// now fill in the M and operator etc.
	O = sq.getOper();
	mvalue = sq.getM();
	division = sq.getDivision();

	// exp factor
	expFacListLen = sq.getExpFacListLen();
	for(int i=0; i<MAX_EXP_FAC_LIST; i++) expFacList[i] = NULL_POS;
	if (expFacListLen>0) {
		for(int i=0; i<expFacListLen; i++) {
			expFacList[i] = sq.getExpFacVal(i);
		}
	}

	// derivative information
	firstDerivPos  = sq.get1stDerivPos();  
	secondDerivPos = sq.get2edDerivPos(); 
	firstDerivDir  = sq.get1stDerivDir();   
	secondDerivDir = sq.get2edDerivDir();  
}

const Basis& Integral::getBasis(const int& pos) const {
	if (pos == BRA1) {
		return bra1;
	}else if (pos == BRA2) {
		return bra2;
	}else if (pos == KET1) {
		return ket1;
	}else {
		crash(pos != KET2,"Wrong of basis passed in getBasis in integral class");
		return ket2;
	}
}

bool Integral::operator==(const Integral& I) const {
	if(O != I.O || bra1 != I.bra1 || bra2 != I.bra2 ||
			ket1 != I.ket1 || ket2 != I.ket2 || 
			mvalue != I.mvalue || division != I.division) { 
		return false;
	}

	// compare deriv infor
	if (firstDerivPos != I.firstDerivPos || secondDerivPos != I.secondDerivPos 
			|| firstDerivDir != I.firstDerivDir || secondDerivDir != I.secondDerivDir) {
		return false;
	}

	// if this integral has different exp factor length 
	// with input one, of course they are not same
	if (expFacListLen != I.expFacListLen) return false;

	// now we need to compare the expotential factors
	if (expFacListLen>0) {
		const int* expList = I.getExpFacList();
		for(int i=0; i<expFacListLen; i++) {
			if (expFacList[i] != expList[i]) return false;
		}
	}
	return true;
}

string Integral::getName() const {
	string name = "I";
	string Oname = getOperStringName(O);
	name = name + "_" + Oname;
	name = name + "_" + bra1.getName();
	if (!bra2.isnull())
		name = name + "_" + bra2.getName();
	if (!ket1.isnull())
		name = name + "_" + ket1.getName();
	if (!ket2.isnull())
		name = name + "_" + ket2.getName();
	if (mvalue > 0)
		name = name + "_" + "M" + lexical_cast<string>(mvalue);
	if (division >= 0)
		name = name + "_" + "C" + lexical_cast<string>(division);

	// now let's add in deriv information
	int jobOrder = 0;
	if (firstDerivPos != NULL_POS) jobOrder++;
	if (secondDerivPos != NULL_POS) jobOrder++;
	for(int i=1; i<=jobOrder; i++) {
		int pos = firstDerivPos; 
		if (i == 2) pos = secondDerivPos;
		if (pos > 0) {
			name = name + "_d";
			if (pos == BRA1) {
				name = name + "a";
			}else if (pos == BRA2) {
				name = name + "b";
			}else if (pos == KET1) {
				name = name + "c";
			}else if (pos == KET2) {
				name = name + "d";
			}else{
				crash(true,"in the getName of integral class, when pos is >0 but value is invalid?");
			}
		}
		int dir = firstDerivDir;
		if (i == 2) dir = secondDerivDir;
		if (dir > 0) {
			if (dir == DERIV_X) {
				name = name + "x";
			}else if (dir == DERIV_Y) {
				name = name + "y";
			}else if (dir == DERIV_Z) {
				name = name + "z";
			}else{
				crash(true,"in the getName of integral class, when dir is > 0 but value is invalid?");
			}
		}
	}

	// expotential factor
	if (expFacListLen>0) {
		for(int i=0; i<expFacListLen; i++) {
			if (i == 0) {
				name = name + "_";
			}
			int pos = expFacList[i]; 
			if (pos > 0) {
				if (pos == BRA1) {
					name = name + "a";
				}else if (pos == BRA2) {
					name = name + "b";
				}else if (pos == KET1) {
					name = name + "c";
				}else if (pos == KET2) {
					name = name + "d";
				}
			}
		}
	}
	return name;
}

/*
bool Integral::isSTypeIntegral() const
{
	const Basis& sbra1 = getBasis(BRA1); 
	const Basis& sbra2 = getBasis(BRA2); 
	const Basis& sket1 = getBasis(KET1); 
	const Basis& sket2 = getBasis(KET2); 
	int L = 0;
	L = sbra1.getL();
	int nBody = getOperOrder(O);
	if (! sbra2.isnull() && nBody>=2) L += sbra2.getL();
	if (! sket1.isnull() && nBody>=3) L += sket1.getL();
	if (! sket2.isnull() && nBody>=4) L += sket2.getL();
	if (L == 0) return true;
	return false;
}
*/

string Integral::formVarName(const int& rrType) const
{
	string varName = getName();

	// do we need to add in the modifier of _vrr?
	// all of VRR results (local and module results comes with _vrr)
	if (isValidVRRJob(rrType)) {
		varName = varName + "_vrr";
	}

	return varName;
}

int Integral::getIndex() const {

	// get total number of basis set in this shell - we calculate directly
	// also together with the local index of the basis set
	// they are incremented into the final index value
	
	// bra1
	int L = bra1.getL();
	Shell s(L);
	int n1 = s.getBasisSetNumber();
	int a1 = bra1.getLocalIndex();
	int index  = a1;

	// bra2
	int n2 = 0;
	int a2 = 0;
	if (! bra2.isnull()) {
		int L = bra2.getL();
		Shell s(L);
		n2 = s.getBasisSetNumber();
		a2 = bra2.getLocalIndex();
		index += a2*n1;
	}

	// ket1
	int n3 = 0;
	int a3 = 0;
	if (! ket1.isnull()) {
		int L = ket1.getL();
		Shell s(L);
		n3 = s.getBasisSetNumber();
		a3 = ket1.getLocalIndex();
		index += a3*n1*n2;
	}

	// ket2
	// we difinitely do not need the dimension for ket2
	int a4 = 0;
	if (! ket2.isnull()) {
		a4 = ket2.getLocalIndex();
		index += a4*n1*n2*n3;
	}

	return index;
}

void Integral::addExpFac(int pos) 
{
	// first step, we need double check
	crash(expFacListLen>=MAX_EXP_FAC_LIST || expFacListLen<0, 
			"invalid exp fac length in integral class, already >=MAX_EXP_FAC_LIST or < 0");
	if (pos != BRA1 && pos != BRA2 && pos != KET1 && pos != KET2) {
		crash(true, "invalid position value pass in addExpFac in integral class");
	}

	// now add in value
	expFacList[expFacListLen] = pos; 
	expFacListLen++;

	// we need to re-shuffle the value so that to make it 
	// in order
	int list1[MAX_EXP_FAC_LIST];  
	for(int i=0; i<MAX_EXP_FAC_LIST; i++) list1[i] = expFacList[i];
	std::sort(list1, list1+MAX_EXP_FAC_LIST);

	// reset the expFacList
	for(int i=0; i<MAX_EXP_FAC_LIST; i++) expFacList[i] = NULL_POS;

	// now let's copy it back
	// just omit the NULL values
	int j=0;
	for(int i=0; i<MAX_EXP_FAC_LIST; i++) {
		if (list1[i] == NULL_POS) continue;
		expFacList[j] = list1[i];
		j++;
	}
}

