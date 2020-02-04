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

#ifndef _GAUSS_MATRICES_H_
#define _GAUSS_MATRICES_H_

#include "Matrix.h"
#include "../Integration/IntegrationScheme.h"

namespace FluxSol
{


	class GaussFullMatrices:
		public MatrixIndex
	{

		protected:
			GaussIntegrationScheme intsch;	//CONST???
			vector< Matrix<double> > vm;
		public:
			GaussFullMatrices(){}
			GaussFullMatrices(const int &r, const int &c, const GaussIntegrationScheme &is)
			:MatrixIndex(r, c),
			intsch(is) //Has constant members
			{
				intsch = is;
				for (int i = 0; i < is.NumPoints(); i++)
					vm.push_back(Matrix<double>(r, c));
					//cout << "is numpoints"<<is.NumPoints()<<endl;
			}
			inline const Matrix<double> Mat(const unsigned int &gausspoint)const{ return vm[gausspoint]; };

			//NOT CONST since is modified by these
			inline Matrix<double> & operator[](const int &i){ return vm[i]; }

			inline const std::string outstr()const;	//TO MODIFY: MAKE VIRTUAL
			const GaussIntegrationScheme & IntScheme()const { return intsch; }

			const GaussFullMatrices operator *(const Matrix<double> &m2) const;
			//NOT CONSt, TO MODIFY

			inline GaussFullMatrices& operator=(const GaussFullMatrices &m)
			{
				this->rows = m.Rows();
				this->cols = m.Cols();
				GaussFullMatrices m2(m);
				for (int s = 0; s<m2.vm.size();s++)
				this ->vm.push_back(m2.Mat(s));
				this->intsch = m2.intsch;
				return *this;
			}

			const double & cval(const int &g, const int &r, const int &c)const { return this->vm[g].cval(r,c); }


	};
	//Gauss  matrix evaluated in every gauss point
	//Elemental matrix with
	class GaussMatrices:
		public SparseMatrixIndex
	{
		GaussIntegrationScheme intsch;	//CONST???
	protected:

		std::vector<SparseMatrixBulk> vm;
		//const unsigned int valnumber;

	public:
		//MUST INITIALIZE DIMENSION
		inline GaussMatrices(){};
		inline GaussMatrices(const int &rows, const int &cols, const GaussIntegrationScheme &intsch);
		inline GaussMatrices(const int &rows, const int &cols, const int &numvals,const GaussIntegrationScheme &intsch);
		inline const std::vector< Matrix<double> > Mat()const;
		inline void SetPair(const unsigned int &p, const unsigned int &row, const unsigned int &col)
		{ this->pos[p].Set(row,col); }
		inline const Matrix<double> Mat(const unsigned int &gausspoint)const;

		const std::string outstr()const;
		inline const std::string outstr(const int &gausspoint)const;




		inline SparseMatrixBulk& operator[](const int &i){ return vm[i]; }


		//MATH OPERATORS
		inline const GaussMatrices operator *(const GaussMatrices &)const;
		inline const GaussMatrices operator +(const GaussMatrices &)const;

	};


	class ShapeFunctionMatrices :	//SEVERAL GAUSS VALUES
		//public FEGPValuedMatrices,
		public GaussMatrices
	{
		//std::vector <ShapeFunction> shapefn;		//vector of linear functions, depends on element type
		//DOES NOT MAKE SENSE TO MAINTAIN THIS IN MEMORY
		//ShapeFunctionGroup shapefngrp;				//reference??

	public:
		ShapeFunctionMatrices();
		//ShapeFunctionMatrices(/*const Element &e*/);
		//The next statement can be done in a function
		//Matrix <double> rst(const Vec3D &);	//evaluate at rst value
		Matrix <double> GaussPoint(const int &);	//evaluate at gauss point index

	};


	//INLINE FUNCTIONS

	//-------------
	// CONSTRUCTORS
	//-------------

	GaussMatrices::GaussMatrices(const int &rows, const int &cols, const GaussIntegrationScheme &intsch)
	{

		for (int p = 0; p < intsch.NumPoints(); p++)
		{
			this->vm.push_back(SparseMatrix(rows, cols));

		}


	}

	inline GaussMatrices::GaussMatrices(const int &rows, const int &cols, const int &numvals, const GaussIntegrationScheme &intsch):
		SparseMatrixIndex(rows,cols)
	{
		for (int p = 0; p < intsch.NumPoints(); p++)
			this->vm.push_back(SparseMatrixBulk(numvals));

		for (int n=0; n < numvals; n++)
			this->pos.push_back(OrderedPair());
	}

	const std::vector< Matrix<double> > GaussMatrices::Mat() const
	{
		vector< Matrix<double> > vm;

		Matrix<double> Mat;


		return vm;
	}



	const Matrix<double> GaussMatrices::Mat(const unsigned int &gausspoint)const
	{
		Matrix<double> m(this->rows, this->cols);

		for (unsigned pos = 0; pos < this->pos.size(); pos++)
			m[this->Pair(pos).Row()][this->Pair(pos).Col()] = vm[gausspoint].Val(pos);

		return m;
	}


	inline const std::string GaussMatrices::outstr()const
	{
		std::string s;
		Matrix<double> m;

		for (int gp = 0; gp < vm.size(); gp++)
		{
			s += "Gauss Point \n";
			//s += this->Mat(gp).outstr();
			m = this->Mat(gp);
			s += m.outstr();
		}

		return s;
	}

	//TO MODIFY: MAKE VIRTUAL
	inline const std::string GaussFullMatrices::outstr()const
	{
		std::string s;
		Matrix<double> m;

		for (int gp = 0; gp < vm.size(); gp++)
		{
			s += "Gauss Point \n";
			//s += this->Mat(gp).outstr();
			m = this->Mat(gp);
			s += m.outstr();
		}

		return s;
	}






	inline const GaussFullMatrices GaussFullMatrices::operator *(const Matrix<double> &m2)const
	{
		const GaussFullMatrices m1 = *this;

		GaussFullMatrices m(m1.Rows(), m2.Cols(), m1.IntScheme());

        cout << "Num Points "<<m1.IntScheme().NumPoints()<<endl;
		for (int g = 0; g<m1.IntScheme().NumPoints(); g++)
		if (m1.Cols() == m2.Rows())
		{
			for (int f = 0; f < m1.Rows(); f++)
			{
				for (int cext = 0; cext < m2.Cols(); cext++)
				{

                    //cout << "g f crxt "<< g << " " << f << " "<< cext << " "<<endl;
					double val = 0.0;
					for (int c = 0; c < m1.Cols(); c++)
						val += m1.cval(g,f,c) * m2.cval(c,cext);

					m[g][f][cext] = val;
				}

			}
		}


		return m;

	}

};//FluxSol

#endif
