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
 * \file    jdigestion.h
 * \brief   describing the unrolling of j digestion
 * \author  Fenglai Liu
 */
#ifndef JDIGESTION_H
#define JDIGESTION_H
#include "general.h"
#include "shellquartet.h"
using namespace shellquartet; 

namespace jdigestion {


	class JDigestion {

		private:

			ShellQuartet sq;

		public:

			///
			/// constructor
			///
			JDigestion(const ShellQuartet& sq0):sq(sq0) { };

			///
			/// destructor
			///
			~JDigestion() = default;

            ///
            /// initialize the string vector which holds the result
            /// we will put the result block name as initialization
            ///
            void initResultVectionString(bool isBraSide, vector<string>& result) const;

			///
			/// for the unrolling, get either density matrix value name;
			/// or result name
			///
			string get2DBlockName(bool isDenMtrxBlock,
					const Basis& b1, const Basis& b2, bool isBraSide) const;

			///
			/// unrolling the J digestion
			///
			void unrollingJDigestion() const;
	};
}

#endif
