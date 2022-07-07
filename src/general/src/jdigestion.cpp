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
#include "jdigestion.h"
using boost::lexical_cast;
using namespace basisutil;
using namespace basis;
using namespace shell;
using namespace integral;
using namespace jdigestion;

string JDigestion::get2DBlockName(bool isDenMtrxBlock, const Basis &b1,
                                  const Basis &b2, bool isBraSide) const {

    // prefix
    string name = "result_";
    if (isDenMtrxBlock) name = "dm_";

    // add in basis set name
    name += b1.getName() + "_" + b2.getName();

    // add in pos infor
    if (isBraSide) {
        name += "_bra";
    } else {
        name += "_ket";
    }

    // return
    return name;
}

void JDigestion::initResultVectionString(bool isBra, vector<string>& result) const {

    // get the shell data
    BasisUtil bu;
    const Shell& bra1 = sq.getShell(BRA1);
    const Shell& bra2 = sq.getShell(BRA2);
    const Shell& ket1 = sq.getShell(KET1);
    const Shell& ket2 = sq.getShell(KET2);
    int nBra1 = bra1.getBasisSetNumber();
    int nBra2 = bra2.getBasisSetNumber();
    int nKet1 = ket1.getBasisSetNumber();
    int nKet2 = ket2.getBasisSetNumber();

    // set the row and column
    int nRow = nKet1;
    int nCol = nKet2;
    int rowL = ket1.getL();
    int colL = ket2.getL();
    if (isBra) {
        nRow = nBra1;
        nCol = nBra2;
        rowL = bra1.getL();
        colL = bra2.getL();
    }

    // set the length
    result.assign(nRow*nCol," ");

    // now generate the initialization
    for (int j=0; j<nCol; j++) {

        // get the basis name
        int l2,m2,n2;
        bu.getLMNVal(colL,j,l2,m2,n2);
        Basis k1(l2,m2,n2);

        // loop over the row basis set
        for (int i=0; i<nRow; i++) {

            // get the basis name
            int l1,m1,n1;
            bu.getLMNVal(rowL,i,l1,m1,n1);
            Basis b1(l1,m1,n1);

            // get the result block name
            string name = get2DBlockName(false, b1, k1, isBra) + " = ";

            // get the index
            int index = i+j*nRow;
            result[index] = name;
        }
    }
}

void JDigestion::unrollingJDigestion() const {

    // set up vectors to store the results
    vector<string> j12;
    vector<string> j34;
    initResultVectionString(true,  j12);
    initResultVectionString(false, j34);

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

            // let's get the ket side index
            int ketIndex = k+l*nKet1;

            // density matrix value
            string p34 = get2DBlockName(true, k1, k2, false);

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

                    // bra side index
                    int braIndex = i+j*nBra1;

                    // now we can get the integral
                    Integral I(b1,b2,k1,k2,sq.getOper());
                    string iname = I.getName();

                    // p12
                    string p12 = get2DBlockName(true, b1, b2, true);

                    // now add the line
                    string j12Line = iname + "*" + p34 + " + ";
                    j12[braIndex] += j12Line;

                    // ket side
                    string j34Line = iname + "*" + p12 + " + ";
                    j34[ketIndex] += j34Line;
                }
            }
        }
    }

    // finally let's print out the result
    for(int i=0; i<j12.size(); i++) {
        printf("%-s\n", j12[i].c_str());
    }
    printf("\n\n");
    for(int i=0; i<j34.size(); i++) {
        printf("%-s\n", j34[i].c_str());
    }
    printf("\n\n");
}

