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
    string name = "result_";
    if (isDenMtrxBlock) name = "dm_";

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

int KDigestion::get2DMatrixIndex(const int& i, const int& j, const int& pos) const {

    // get the shell
    const Shell& bra1 = sq.getShell(BRA1);
    const Shell& bra2 = sq.getShell(BRA2);
    int nBra1 = bra1.getBasisSetNumber();
    int nBra2 = bra2.getBasisSetNumber();

    // get the row dimension
    int nRowBas = nBra1;
    if (pos == K23 || pos == K24) nRowBas = nBra2;

    // result
    int index = i+j*nRowBas;
    return index;
}

void KDigestion::initResultVectionString(const int& pos, vector<string>& result) const {

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
    int nRow = nBra1;
    int nCol = nKet1;
    int rowL = bra1.getL();
    int colL = ket1.getL();
    if (pos == K23) {
        nRow = nBra2;
        rowL = bra2.getL();
    } else if (pos == K14) {
        nCol = nKet2;
        colL = ket2.getL();
    } else if (pos == K24) {
        nRow = nBra2;
        nCol = nKet2;
        rowL = bra2.getL();
        colL = ket2.getL();
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
            string name = get2DBlockName(false, b1, k1, pos) + " = ";

            // get the index
            int index = get2DMatrixIndex(i,j,pos);
            result[index] = name;
        }
    }
}

void KDigestion::unrollingKDigestion() const {

    // set up vectors to store the results
    vector<string> k13;
    vector<string> k23;
    vector<string> k14;
    vector<string> k24;
    initResultVectionString(K13, k13);
    initResultVectionString(K23, k23);
    initResultVectionString(K14, k14);
    initResultVectionString(K24, k24);

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

                    // compute the index for the result
                    int k13_index = get2DMatrixIndex(i,k,K13);
                    int k14_index = get2DMatrixIndex(i,l,K14);
                    int k23_index = get2DMatrixIndex(j,k,K13);
                    int k24_index = get2DMatrixIndex(j,l,K14);

                    // now we can get the integral
                    Integral I(b1,b2,k1,k2,sq.getOper());
                    string iname = I.getName();

                    // k13 - p24
                    string p24 = get2DBlockName(true, b2, k2, K24);

                    // k23 - p14
                    string p14 = get2DBlockName(true, b1, k2, K14);

                    // k14 - p23
                    string p23 = get2DBlockName(true, b2, k1, K23);

                    // k24 - p13
                    string p13 = get2DBlockName(true, b1, k1, K13);

                    // now let's print
                    string line1 = iname + "*" + p24 + " + ";
                    k13[k13_index] += line1;
                    string line2 = iname + "*" + p14 + " + ";
                    k23[k23_index] += line2;
                    string line3 = iname + "*" + p23 + " + ";
                    k14[k14_index] += line3;
                    string line4 = iname + "*" + p13 + " + ";
                    k24[k24_index] += line4;
                }
            }
        }
    }

    // finally let's print out the result
    for(int i=0; i<k13.size(); i++) {
        printf("%-s\n", k13[i].c_str());
    }
    printf("\n\n");
    for(int i=0; i<k23.size(); i++) {
        printf("%-s\n", k23[i].c_str());
    }
    printf("\n\n");
    for(int i=0; i<k14.size(); i++) {
        printf("%-s\n", k14[i].c_str());
    }
    printf("\n\n");
    for(int i=0; i<k24.size(); i++) {
        printf("%-s\n", k24[i].c_str());
    }
}

