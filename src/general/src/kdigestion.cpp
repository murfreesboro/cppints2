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
#include "boost/lexical_cast.hpp"
#include "basisutil.h"
#include "basis.h"
#include "shell.h"
#include "integral.h"
#include "kdigestion.h"
using boost::lexical_cast;
using namespace basisutil;
using namespace basis;
using namespace shell;
using namespace integral;
using namespace kdigestion;

string KDigestion::get2DBlockName(bool isDenMtrxBlock, const Basis &b1,
                                  const Basis &b2, const int& pos) const {

    // prefix
    string name = "result";
    if (isDenMtrxBlock) name = "dm";

    // add in basis set name
    name += b1.getName() + "_" + b2.getName();

    // add in pos infor
    if (pos == K13) name += "_k13";
    if (pos == K14) name += "_k14";
    if (pos == K23) name += "_k23";
    if (pos == K24) name += "_k24";

    // return
    return name;
}

void KDigestion::unrollingKDigestion() const {

    // get the shells
    BasisUtil bu;
    const Shell& bra1 = sq.getShell(BRA1);
    const Shell& bra2 = sq.getShell(BRA2);
    const Shell& ket1 = sq.getShell(KET1);
    const Shell& ket2 = sq.getShell(KET2);
    int nBra1 = bra1.getBasisSetNumber();
    int nBra2 = bra2.getBasisSetNumber();
    int nKet1 = ket1.getBasisSetNumber();
    int nKet2 = ket2.getBasisSetNumber();

    // let's loop
    for(int l=0; l<nKet2; l++) {

        // get the basis set
        int l4,m4,n4;
        bu.getLMNVal(ket2.getL(),l,l4,m4,n4);
        Basis k2(l4,m4,n4);

        for(int k=0; k<nKet1; k++) {

            // get the basis set
            int l3,m3,n3;
            bu.getLMNVal(ket1.getL(),k,l3,m3,n3);
            Basis k1(l3,m3,n3);

            for (int j=0; j<nBra2; j++) {

                // get the basis set
                int l2,m2,n2;
                bu.getLMNVal(bra2.getL(),j,l2,m2,n2);
                Basis b2(l2,m2,n2);

                for(int i=0; i<nBra1; i++) {

                    // get the basis set
                    int l1,m1,n1;
                    bu.getLMNVal(bra1.getL(),i,l1,m1,n1);
                    Basis b1(l1,m1,n1);

                    // now we can get the integral
                    Integral I(b1,b2,k1,k2,sq.getOper());
                    string iname = I.getName();

                    // k13 - p24
                    string k13 = get2DBlockName(false, b1, k1, K13);
                    string p24 = get2DBlockName(true, b2, k2, K24);

                    // k23 - p14
                    string k23 = get2DBlockName(false, b2, k1, K23);
                    string p14 = get2DBlockName(true, b1, k2, K14);

                    // k14 - p23
                    string k14 = get2DBlockName(false, b1, k2, K14);
                    string p23 = get2DBlockName(true, b2, k1, K23);

                    // k24 - p13
                    string k24 = get2DBlockName(false, b2, k2, K24);
                    string p13 = get2DBlockName(true, b1, k1, K13);

                    // now let's print
                    string line1 = "Double " + k13 + " += " + iname + "*" + p24 + ";";
                    printf("%-s\n", line1.c_str());
                    string line2 = "Double " + k23 + " += " + iname + "*" + p14 + ";";
                    printf("%-s\n", line2.c_str());
                    string line3 = "Double " + k14 + " += " + iname + "*" + p23 + ";";
                    printf("%-s\n", line3.c_str());
                    string line4 = "Double " + k24 + " += " + iname + "*" + p13 + ";";
                    printf("%-s\n", line4.c_str());
                }
            }
        }
    }
}

