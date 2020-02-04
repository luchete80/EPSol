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

#ifndef _MATERIALMATRIX_H_
#define _MATERIALMATRIX_H_


#include "../../Type/Matrix.h"

namespace FluxSol{

	typename <int dim>
	class MaterialMatrix:public Matrix<double>
	{
	

	protected:
		virtual void Create();

	public:

		MaterialMatrix(){ Create(); }

	};

	class LinearElasticMaterialMatrix
		:public MaterialMatrix
	{

		public:

			//Constructors
			LinearElasticMaterialMatrix() : MaterialMatrix(){}
			void Create();	//virtual function

	};

	class PlainStrainLinearElasticMaterialMatrix
	:public MaterialMatrix
	{

		public:

			//Constructors
			LinearElasticMaterialMatrix() : MaterialMatrix(){}
			void Create();	//virtual function

	};


}//FluxSol

#endif
