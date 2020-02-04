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

#ifndef _FEMATRIX_H_
#define _FEMATRIX_H_

#include "Matrix/Matrix.h"

#include "../Mesh/Element.h"
#include "../ShapeFunction/ShapeFunction.h"
#include "Integration\IntegrationScheme.h" //BUT MUST NOT INCLUDE FEINTEGRATION SINCE IT INCLUDES THIS FILE

//#include "../Integration/IntegrationScheme.h"

#include <vector>

using namespace std;


namespace FluxSol
{



	//ThESE aRE not REALLY matriceS, but returns them


	//THIS CLASS INCLUDES EVERY FIELD, IS THE GLOBAL MATRIX
	//DONT MAKE ANY SENSE TO KEEP ALL VALUE SINCE SEVERAL ARE ZEROS
	class FESparseMatrix
		:public SparseMatrix
	{

	protected:

	public:
		FESparseMatrix(){};


	};

	template <int dim>
	class FEElemMatrix:
		public FESparseMatrix
	{
	protected:

		//friend class fvSolver;
		const Element <dim> &elref;	//Constant reference to element
	public:
		FEElemMatrix(const Element <dim> &ref) :elref(ref){}
		FEElemMatrix();
	};

	//References to an element
	template <int dim>
	class FEElemSparseMatrix :
		public FESparseMatrix,
		public FEElemMatrix<dim>
	{

		//- Const reference to GeometricField<Type, fvPatchField, volMesh>
		//  Converted into a non-const reference at the point of solution.
		//const GeometricField<Type, fvPatchField, volMesh>& psi_;
	protected:


		FEElemSparseMatrix():FEElemMatrix<dim> (Element<dim>()){}
		FEElemSparseMatrix(const Element <dim> &ref) :FEElemMatrix<dim>(ref){ }

	public:

	};

	//THIS CLASS INCLUDES EVERY FIELD, IS THE GLOBAL MATRIX
	class FESqSparseSymmMatrix
		:public FESparseMatrix
	{

	protected:

	public:


	};


	//THIS CLASS INCLUDES EVERY FIELD, IS THE GLOBAL MATRIX
	class FEGlobalMatrix
		:public FESqSparseSymmMatrix
	{

	protected:

	public:


	};


	//A finite element matrix evaluated in every gauss point
	//Elemental matrix with
	class FEGaussElemMatrices
	{
		GaussIntegrationScheme intsch;	//CONST???
	protected:
		std::vector<FESparseMatrix> vm;

	public:
		inline FEGaussElemMatrices(){};
		inline FEGaussElemMatrices(const int &rows, const int &cols, const GaussIntegrationScheme &intsch);

	};



	//Single Element Linear StiffNessMatrix
	//SINGLE VALUE MATRIX
	template <int dim>
	class LinearStiffnessFEElemMatrix:
		public FESqSparseSymmMatrix,
		public FEElemMatrix<dim>
	{


	public:
		inline LinearStiffnessFEElemMatrix(const Element <dim> &e);
	};




	//LinearStiffnessFEElemMatrix::LinearStiffnessFEElemMatrix(const Element &e)
	//	:FEElemMatrix(e)//Initializes reference
	//{
	//

	//	ShapeFunctionMatrices H(e);


	//
	//	//Derives ShapeFunction Matrix
	//	//Diff(H,Coords(x0)) where x0 is reference system
	//
	//}


	//INLINE FUNCTIONS

	FEGaussElemMatrices::FEGaussElemMatrices(const int &rows, const int &cols, const GaussIntegrationScheme &intsch)
	{


	}


};




#endif
