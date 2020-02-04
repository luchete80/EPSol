/************************************************************************

Copyright 2012-2013 Luciano Buglioni

Contact: luciano.buglioni@gmail.com

This file is a part of FluxSol

FluxSol is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

Free CFD is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

For a copy of the GNU General Public License,
see <http://www.gnu.org/licenses/>.

*************************************************************************/

#ifndef _FUNCTION_H_
#define _FUNCTION_H_

#include <vector>
using namespace std;

#include "Matrix/CoeffMatrix.h"

//Function of different type return
namespace FluxSol
{

	template <typename T>
	class Function {

	protected:

		//Linear terms, always are present;
		double kterm;		//Constant term
		std::vector <double> linearcoeff;	//Six terms, r,s,t and cross



	public:
		Function(){};


	};

	//Proportional to variables function, i.e., these are multiplying
	template <typename T>
	class PropFunction :
		public Function <double>{

	protected:

		//Linear terms, always are present;
		T kterm;		//Constant term
		std::vector <T> linearcoeff;	//Six terms, r,s,t and cross

		CoeffMatrix coeff;				//GENERAL COEFFICIENT MATRIX
		int order_; // TO MODIFY: CONST


	public:
		PropFunction():order_(0){};
		PropFunction(const int &order):coeff(order,order,order),order_(order){};
		inline const std::vector<PropFunction> diff() const;
		const int& Order()const{ return order_; }

	};




	//INLINE FUNCTIONS
	template <typename T>
	const std::vector< PropFunction<T> > PropFunction<T>::diff() const
	{
		PropFunction<T> fr, fs, ft;

		//If exponent>=0, derivate
		for (int var = 0; var < order_; var++)
		{
			for (int r = 1; r<order_; r++)	fr.coeff[r - 1][var][var] += r*this->coeff.Val(r,var,var);
			for (int s = 1; s<order_; s++)	fs.coeff[var][s - 1][var] += s*this->coeff.Val(var,s,var);
			for (int t = 1; t<order_; t++)	ft.coeff[var][var][t - 1] += t*this->coeff.Val(var,var,t);
		}

		vector< PropFunction<T> > vf;
		vf.push_back(fr, fs, ft);


		return vf;
	}


}//FluxSol
#endif //Function
