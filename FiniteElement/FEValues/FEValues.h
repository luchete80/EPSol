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

#ifndef _FE_VALUES_H
#define _FE_VALUES_H

#include <iostream>

#include "Element.h"
#include "../FEMatrix/FeMatrix.h"
#include "../Type/Vec3d.h"
#include "FEGrid.h"
#include "../../Matrix/GaussMatrices.h"

#include "Type/Operations.h"

#include <vector>

using namespace std;

namespace FluxSol
{

	//TO MODIFY int dim, int number
	template <int dim>
	class FEValues{

	protected:
		bool modify;

		//TO MODIFY, TO MAKE A TEMPLATE

		const Element<dim> &elref;


		const FeGrid<dim> &grid;

		//ITERATOR LIKe in dealii???

		//TEMPORARY DATA
		GaussMatrices temp;
		GaussFullMatrices shape_value_matrices;
		GaussFullMatrices jacobian;
		GaussFullMatrices shape_grad_matrices;			//Matrices with linear arrangements to do Bt x C x B
		GaussFullMatrices shape_localgrad_matrices;

		GaussFullMatrices shape_grad_components;		//Without zeroes
		GaussFullMatrices shape_value_components;		//Without zeroes, 0 x nodenumber matrix


		//vector<Vector<double>> ShapeFnVals
		Matrix<double> xyz;

	public:

		//Options for initialize
		FEValues():elref(Element<dim>()){
			this->SetJacobian;
		}

		//TO MODIFY: ONLY NEED ELEMENT
		FEValues(const Element<dim> &e):elref(e),grid(FeGrid<dim>())
		{

		}

		FEValues(const Element<dim> &e, const FeGrid<dim> &g) :elref(e), grid(g){

            cout << "Setting Shape values .."<<endl;
			this->Set_shape_value_matrix();
            cout << "Setting Shape values .."<<endl;
			this->SetJacobian();
            cout << "Setting Shape values .."<<endl;
			this->Set_shape_grad_matrix();

		}
		//OR VECTOR
		//MUST BE Vector EPSol class
		// TO MODIFY
		inline void SetJacobian();
		inline void SetShapeValue();
		inline void SetShapeGrad();

		inline const GaussFullMatrices & Jacobian()const;
		//TO MODIFY, TENSOR
		const vector<double> & shape_grad_component(const int &f, const int &gausspoint)const;
		const Matrix<double> & shape_grad_comps(const int &g) const{return this->shape_grad_components.Mat(g);}
		const GaussFullMatrices & shape_grad_comps() const{return this->shape_grad_components;}
		const Matrix<double> & shape_value_matrix(const int &gausspoint)const;
		inline const GaussFullMatrices & shape_value_matrix()const{return this->shape_value_matrices;}
		inline const double & shape_grad_component  (const int &fn, const int &gaussp, const int &comp) const{return this->shape_grad_components[gaussp][comp][fn];};	//Like in dealii
		inline const double & shape_value_component (const int &fn, const int &gaussp) const{return this->shape_value_components[gaussp][0][fn];};	//Like in dealii
		//TO MODIFY TO shape_grad_matrix
		inline const GaussFullMatrices  & shape_grad_matrix()const;
		inline const GaussFullMatrices  & shape_localgrad_matrix()const;

		inline const Matrix<double> MatMatrix()const{cout << "returning"<<endl;return elref.MatMatrix();}	//TO MODIFY, SAVE MATRIX IN CLASS

		//Set value gauss Matrices
		inline void virtual Set_shape_value_matrix();
		inline void virtual Set_shape_grad_matrix();
		const unsigned int DoF_Number() const {return elref.DoFpNode();}


	};//FEValues



	template <int dim>
	inline void FEValues<dim>::SetJacobian()
	{
		GaussIntegrationScheme gi(elref.GaussOrder(), dim);
		cout << "GaussOrder "<<elref.GaussOrder()<<endl;
		//GaussIntegrationScheme gi();
		//Or elref.DIM()//
		GaussFullMatrices ret(dim, dim, gi);
		cout << "GaussOrder gi"<<gi.NumPoints()<<endl;

		ShapeFunctionGroup shfngr = elref.CreateShapeFunctionGroup();

		GaussFullMatrices temph(dim, shfngr.Size(), gi);
		//vector<ShapeFunction> vsh;
		//Evaluating in Gauss Points
		//double temp;
		for (int g = 0; g < gi.NumPoints(); g++)
		{
			for (int h = 0; h < shfngr.Size(); h++)
			{
				const vector<ShapeFunction> vsh = shfngr.ShapeFn(h).diff();
					for (int comp = 0; comp < dim; comp++)
					{
						//temp=vsh[comp].Val(gi[g]);
						temph[g][comp][h] = vsh[comp].Val(gi[g]);
						//temph[g][comp][h] = 1.;
					}

			}

			//	//can be done outside for
			//	//can be done outside temph*this->grid.XYZ(elref);

			//const


		}
		//[dh1dr dh2dr ..]   x   [X1 Y1]
		//[dh1ds dh2ds ..]       [X2 Y2]
		//                       [.. ..]
		this->shape_localgrad_matrices = temph;
		//cout << "XYZ"<<this->grid.XYZ(elref).outstr();
		//cout << "temph"<<temph.outstr();
		GaussIntegrationScheme gitemp=temph.IntScheme();
		//cout << "temph int scheme points"<<temph.IntScheme().NumPoints()<<endl;
		this->jacobian = temph*this->grid.XYZ(elref);
		//for (int g = 0; g < gi.NumPoints(); g++)
        //    this->jacobian.Mat(g).Set_Determinant();

	}


	//void Jacobian();
	template <int dim>
	inline const GaussFullMatrices & FEValues<dim>::Jacobian()const
	{

		return jacobian;

	}

	//EXPLICIT
	template <int dim>
	//TO MODIFY TO GaussMatrices
	//GaussMatrices  FEValues<dim>:: shape_grad_matrix()const
	//ASSUME JACOBIAN HAVE BEEN CALCULATED
	//[]
	inline void FEValues<dim>::Set_shape_grad_matrix()
	{

		GaussFullMatrices m;
		GaussIntegrationScheme gi(elref.GaussOrder(), dim);

		ShapeFunctionGroup shfngr = elref.CreateShapeFunctionGroup();

		//GaussFullMatrices temph(2, shfngr.Size(), gi);
        cout << "ceating retcomp gradmatrix"<<endl;
        cout << "col size" << shfngr.Size();
		GaussFullMatrices retcomp(dim, shfngr.Size(), gi);

		//Or elref.DIM()//
		//TO MODIFY, DEPENDS ON VARIABLE DIM
		GaussFullMatrices ret(dim, shfngr.Size()*dim, gi);

		//TO MODIFY, MAKE VIRTUAL
		if (dim==2 && elref.Field_Dim()==2)
		{
			for (int g = 0; g < gi.NumPoints(); g++)
			{
				Matrix<double> dHdrst_T = this->shape_localgrad_matrices.Mat(g).Tr();
				Matrix<double> dHdX = dHdrst_T *this->jacobian.Mat(g).inv();


				for (int h = 0; h < shfngr.Size(); h++)
				{
					//TO MODIFY
					for (int comp = 0; comp < dim; comp++)
					{
						ret    [g][comp][dim * h + comp] = dHdX[h][comp];
						retcomp[g][comp][h] = dHdX[h][comp];
					}

				}

				for (int h = 0; h < shfngr.Size(); h++)
                    //for (int cross=0;cross<)
				{
					ret[g][dim][dim * h    ]=dHdX[h][1];
					ret[g][dim][dim * h + 1]=dHdX[h][0];
				}

			}
			shape_grad_matrices = ret;
			shape_grad_components=retcomp;
		}

		//TO MODIFY, MAKE VIRTUAL
		else if (elref.Field_Dim()==1)
		{
			for (int g = 0; g < gi.NumPoints(); g++)
			{
				Matrix<double> dHdrst_T = this->shape_localgrad_matrices.Mat(g).Tr();
				Matrix<double> dHdX = dHdrst_T *this->jacobian.Mat(g).inv();

                //cout << "dhdrst "<< dHdrst_T.outstr()<<endl;
                //cout << "dhdx "<< dHdX.outstr()<<endl;
                //cout << " Jacobian" << this->jacobian.Mat(g).outstr() <<endl;
                //cout << "Inv Jacobian" << this->jacobian.Mat(g).inv().outstr() <<endl;

				for (int h = 0; h < shfngr.Size(); h++)
				{
					//TO MODIFY
					for (int comp = 0; comp < dim; comp++)
						retcomp[g][comp][h] = dHdX[h][comp];


				}

			}
			shape_grad_matrices = retcomp;
			shape_grad_components=retcomp;
		}
		else
        {
            cout << "ERROR: Shape Grad Matrix Not Set"<<endl;
        }



	}//End Set Shape Grad


	template <int dim>
		inline const GaussFullMatrices & FEValues<dim>::shape_grad_matrix() const
	{
			return this->shape_grad_matrices;
	}


	//Set Shape Value Matrices
	template <int dim>
	//TO MODIFY TO GaussMatrices
	//GaussMatrices  FEValues<dim>:: shape_grad_matrix()const
	//ASSUME JACOBIAN HAVE BEEN CALCULATED
	inline void FEValues<dim>::Set_shape_value_matrix()
	{

		GaussFullMatrices m;
        //cout << " Gauss Order" << elref.GaussOrder()<<endl;
		GaussIntegrationScheme gi(elref.GaussOrder(), dim);
		ShapeFunctionGroup shfngr = elref.CreateShapeFunctionGroup();


		//Or elref.DIM()
		//TO MODIFY

		GaussFullMatrices ret(elref.Field_Dim(), shfngr.Size()*elref.Field_Dim(), gi);
		GaussFullMatrices retcomp(1, shfngr.Size(), gi);

        cout << "Loop through points..."<<endl;
		for (int g = 0; g < gi.NumPoints(); g++) //Gauss Point
		{
		    cout << "point "<<g<<endl;
			for (int h = 0; h < shfngr.Size(); h++) //Shape Function Number
			{
			    cout << "h " << endl;
				//TO MODIFY
				retcomp[g][0][h] = shfngr.ShapeFn(h).Val(gi[g]);
				for (int comp = 0; comp < elref.Field_Dim(); comp++)
					ret[g][comp][(h * elref.Field_Dim() ) + comp] = shfngr.ShapeFn(h).Val(gi[g]);
			}

		}

        cout << "assign matrix"<<endl;
		shape_value_matrices = ret;
        cout << "assign matrix 2"<<endl;
		shape_value_components=retcomp;
		cout << "values assigned "<<endl;
	}

	template <int dim>
	inline const GaussFullMatrices & FEValues<dim>::shape_localgrad_matrix()const
	{
		return this->shape_localgrad_matrices;

	}
}

#endif
