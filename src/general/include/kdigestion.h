/**
 *
 * CPPINTS: A C++ Program to Generate Analytical Integrals Based on Gaussian
 * Form Primitive Functions
 *
 * Copyright (C) 2015 The State University of New York at Buffalo
 * This softare uses the MIT license as below:
 *
 *	Permission is hereby granted, free of charge, to any person obtaining 
 *	a copy of this software and associated documentation files (the "Software"), 
 *	to deal in the Software without restriction, including without limitation 
 *	the rights to use, copy, modify, merge, publish, distribute, sublicense, 
 *	and/or sell copies of the Software, and to permit persons to whom the Software 
 *	is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all 
 * copies or substantial portions of the Software.
 *						    
 *	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
 *	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
 *	PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
 *	FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
 *	ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * \file    kdigestion.h
 * \brief   describing the unrolling of k digestion
 * \author  Fenglai Liu
 */
#ifndef KDIGESTION_H
#define KDIGESTION_H
#include "general.h"
#include "shellquartet.h"
using namespace shellquartet; 

namespace kdigestion {


	class KDigestion {

		private:

			ShellQuartet sq;

		public:

			///
			/// constructor
			///
			KDigestion(const ShellQuartet& sq0):sq(sq0) { };

			///
			/// destructor
			///
			~KDigestion() = default;

            ///
            /// with the input index for the row shell and column shell,
            /// as well as the position; let's compute the index for the
            /// result matrix
            ///
            int get2DMatrixIndex(const int& i, const int& j, const int& pos) const;

            ///
            /// initialize the string vector which holds the result
            /// we will put the result block name as initialization
            ///
            void initResultVectionString(const int& pos, vector<string>& result) const;

			///
			/// for the unrolling, get either density matrix value name;
			/// or result name
			///
			string get2DBlockName(bool isDenMtrxBlock,
					const Basis& b1, const Basis& b2,
					const int& pos) const;

			///
			/// unrolling the K digestion
			///
			void unrollingKDigestion() const;
	};
}

#endif
