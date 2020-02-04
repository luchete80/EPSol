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

#ifndef _SOLVER_H_
#define _SOLVER_H_

namespace FluxSol {


//SQUARE SOLVER
//THIS MAY CONTAIN SYSTEM MATRIX
template <typename number, int dim>
class Solver 
{

	protected:
	
	// Inputs
	number rtol,abstol;
	const int matdim;		//CONST?
	int maxits;
	

	public:
	Solver():matdim(0){}
	Solver (const int &d):
	matdim(d)
	{}
	
	//Virtual destructor
};//Solver

}//FluxSol
#endif