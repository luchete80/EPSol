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
//FINITE ELEMENT INTEGRATION
#ifndef _FEINTEGRATIONSCHEME_H_
#define _FEINTEGRATIONSCHEME_H_

#include "Matrix\Matrix.h"		//DO NOT iNclude FEMATRIX
#include "FiniteElement\FEMatrix\FEMatrix.h"
#include "Integration\IntegrationScheme.h"


//#include "FEField.h"

//INTEGRATION SCHEME DOES NOT SEE MATRIX; MAtrIx referenceS THEM

namespace FluxSol
{


	class FEIntegrationScheme
		:public GaussIntegrationScheme
	{

	public:
		FEIntegrationScheme(){}
		FEIntegrationScheme(const int &i, const int &dim) :GaussIntegrationScheme(i,dim){}

	};



}
#endif
