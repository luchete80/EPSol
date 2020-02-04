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

#ifndef _COEFFMATRIX_H
#define _COEFFMATRIX_H

#include "Matrix.h"

namespace FluxSol{

	//Coefficient Matrix for three variable functions
	class CoeffMatrix{

	protected:
		//Cubic Coefficient
		std::vector< Matrix<double> > v_;
		int dim;

	public:

		//Constructors
		CoeffMatrix(){}
		CoeffMatrix(const int &r_, const int &s_, const int &t_):dim(r_){
			for (int r = 0; r < r_; r++)
				v_.push_back(Matrix<double>(s_, t_));
			for (int r = 0; r < dim;r++)
			for (int s = 0; s < dim; s++)
			for (int t = 0; t < dim; t++)
				this->v_[r][s][t] = 0.0;
		}


		//ACCESS OPERATOR, LIKE IN MATRIX

		//inline const Matrix<double> & operator[](int row)const
		//IT IS NOT CONST SINCE CAN BE USED AS LVALUE
		inline Matrix<double> & operator[](int row)
		{
			return v_[row];
		}

		//USED AS LVALUE
		inline const double & Val(const int&r, const int &s, const int &t)const{ return v_[r].Val(s, t); };


	};

}

#endif
