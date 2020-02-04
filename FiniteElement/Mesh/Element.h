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

#ifndef _ELEMENT_H_
#define _ELEMENT_H_

#include "../../Mesh/Cell.h"

#include "../../Matrix/GaussMatrices.h"

#include "../Mesh/Node.h"
#include "Mat3D.h"

//Virtual functions returning shape functions groups
#include "../ShapeFunction/ShapeFunction.h"

#include "Variables\Variable.h"
#include "Material.h"
//#include "ElType.h"

#include <vector>
#include <fstream>
#include <cstdarg>
#include <iostream>


//TO MODIFY: INCLUDE IN A POLY FUNCTION
namespace FluxSol{

	class Istream;

	template<class PrimitiveType>
	class ElemTraits
	:
		public PrimitiveType
	{

	public:

		// Constructors

			//- Construct from primitive
			ElemTraits(const PrimitiveType& p)
			:
				PrimitiveType(p)
			{}

			//- Construct from Istream
			ElemTraits(Istream& is)
			:
				PrimitiveType(is)
			{}
	};


	//TO MODIFY, INHERIT FROM POLYNOMIAL SPACE
	template <int dim>
	class Element :public _Cell
	{

	protected:
		//std::vector <int> idnode; Is most space demanding
		int num_nodes;	//nodes are not coincident with vertices
//		Mat3D jacobian;	//Or a reference?? In Cell Base class or inherited?? CAN CHANGE
		enum gaussorder{linear=0,quad=1,cubic=2};
		gaussorder gaussor;
		//const int dof;	MEMBER FUNCTION
		//int gaussor;
		unsigned int degree;	//TO MODIFY, INHERIT FROM POLYNOMIAL SPACE
		unsigned int field_dim;	//INHERITED
		int dof;		//TO MODIFY

		Material *mat;

		//returned by function
		//const unsigned int dof_per_node;	//DoF per node

		void AddConnectivity(const std::vector<int> &v);


		//TO SAVE MEMORY NOT ALLOCATED
		// Finite element values
		// GaussFullMatrices shape_value_matrices;
		// GaussFullMatrices jacobian;
		// GaussFullMatrices shape_grad_matrices;			//Matrices with linear arrangements to do Bt x C x B
		// GaussFullMatrices shape_localgrad_matrices;

		// TO MODIFY: possibly not used
		// GaussFullMatrices shape_grad_components;		//Without zeroes
		// GaussFullMatrices shape_value_components;		//Without zeroes, 0 x nodenumber matrix



	public:
		Element()
		{
            gaussor=2;
		}

		Element(const gaussorder &go)
		{
            gaussor=go;
		}

		//Construct from nodes, some of those can be vectex and some others no
		Element (const int &p)	//Polynomial Degree p
		:field_dim(dim)

		{
			this->degree=p;
		}
		Element (const int &p, const gaussorder &go):Element(p){gaussor=go;}
		virtual const unsigned int Field_Dim()const{return field_dim;}//default

		Element(const std::vector<Node> &);
		Element(const std::vector<int> &);
		Element(const int &p, const std::vector<int> &);
		//virtual MakeFromNodes(const std::vector<Node>);

		inline const Matrix<double> XYZ()const;	//return a Matrix with node coordinates to make matrix products
		Mat3D & J(){ return this->jacobian; }
		//TO MODIFY
		const unsigned int GaussOrder()const {
		    unsigned int go=this->gaussor;
		    return go; };

		//inline virtual ShapeFunctionGroup CreateShapeFunctionGroup()const{ return ShapeFunctionGroup(); };
		inline virtual ShapeFunctionGroup CreateShapeFunctionGroup()const;

		inline virtual GaussMatrices H() const;	//Calls
		//inline virtual const GaussFullMatricesHFull()const
		inline virtual GaussMatrices dHdrst() const;	//Calls
		inline virtual const int DOF()const { return 8; }	//TO MODIFY
		inline virtual const int DIM()const { return dim; }

		inline const int & NumNodes()const {return num_nodes; }
		void Set_Nodes(const int numnodes, ...);
		inline void AssignMat(const Material &m){this->mat=&m;};		//TO MODIFY See if Allows inheritance

		const int & NodePos(const int &n)const {return this->id_node[n];}

		//virtual from cell
		std::string Imprimir_Conect(){};

		//inline GaussMatrices GenericJ() const;	//As a 3 x 3 Matrix
		//inline virtual GaussMatrices J() const;

		//inline virtual GaussTensors J() const;

//		inline virtual const FESparseMatrix H(const double &r, const double &s, const double &t) const;	//Calls

		inline const unsigned int & Degree()const{return this->degree;}		//TO MODIFY, INHERIT FROM POLYNOMIAL SPACE
		inline virtual const int DoFPerNode() const {return 0;}

		inline virtual const Matrix<double> MatMatrix()const {cout << "Not virtual"<<endl;};	//Create to save memory


		//----------------------------------- FINITE ELEMENT VALUES ----------------------------------------------
		//TO MODIFY, MAKE PROTECTED
		//
		inline virtual const GaussFullMatrices Shape_grad_matrix() const{};




		//General, not virtual
		inline void Shape_value_matrix();
		//inline const GaussFullMatrices Jacobian() const{return this->jacobian;}
		inline const GaussFullMatrices Jacobian() const;
		inline const GaussFullMatrices Shape_value_matrix()const;

		//Void Functions
		//inline void Set_Jacobian();
	};



	//TYPENAMES

	//PUT ELEMEENT TYPE OR INTEGER
	//typedef QUAD4U   Element<Quad4, DispScheme>;
	//typedef QUAD4UP  Element<Quad4, DispScheme>;
	//typedef QUAD4UPT Element<Quad4, DispScheme>;


	//ASUMING DISPLACEMENT ELEMENT
	//TO DELETE
	class QuadLinearElement
		:public Element<2>
	{


	public:
		QuadLinearElement():Element(){}
		QuadLinearElement(const std::vector<Node> &vn):Element(vn){}

		inline ShapeFunctionGroup CreateShapeFunctionGroup()const;

		inline GaussMatrices H() const{};			//Calls
		inline GaussMatrices dHdrst() const{};	//Calls
		inline const int DOF()const{ return 8; }	//TO MODIFY; THIS DEPENDS ON ELEMENT VARIABLES
		inline const int DIM()const{ return 2; }
		inline const int DoFPerNode() const {return 0;}

	};

	/////////////////////////////////////////// DEFINED ELEMENT TYPES ////////////////////////////////////////////////
	//Scalar Finite Element
	template <int dim>
	class FE_Q:
	public Element<dim>
	{
		public:
		const unsigned int Field_Dim()const {return 1;}	//Scalar dimension
	};

	//Quadrilateral 4 Nodes, plain stress element, integration reduced
	//INCLUDES
	//Q4PSBF: Full Integration
	//Q4PSB
	template <int dim>
	class Q4PS
	:public
	Element<dim>
	{

		public:
		Q4PS(){}
		Q4PS(const std::vector<int> &v)
		:Element<dim>(v)
		{}
		inline const int DOF()const{ return 8; }	//TO MODIFY; THIS DEPENDS ON ELEMENT VARIABLES
		inline const int DIM()const{ return 2; }	//REDUNDANT
		inline const int DoFPerNode() const {return 2;}
		//inline const Matrix<double> MatMatrix const ()
		inline const Matrix<double> MatMatrix () const
		{
			Matrix<double> c(3,3);

			cout << "Creating element Matrix" <<endl;
			double E=this->mat->K(0);
			double nu=this->mat->K(1);

			double ck = E / (1 - nu*nu);
			c[0][0] = c[1][1] = ck;
			c[0][1] = c[1][0] = ck*nu;
			c[2][2] = ck*(1-nu)/2.;

			return c;

			cout << "Virtual created" <<endl;
			cout << c.outstr()<<endl;

		}	//Create to save memory
	};
	//Quadrilateral 4 Nodes, Temperature Element, integration reduced
	template <int dim>
	class Q4TH
	:public
	Element<dim>
	{

	public:
		Q4TH():Element<dim>(){}
        Q4TH(const std::vector<int> &v)
		:Element<dim>(v)
		{
		    //gaussorder d=2;
		     this->gaussor=2;
		}
		Q4TH(const std::vector<int> &v, const unsigned int &go)
		:Element<dim>(v)
		{
		     this->gaussor=go;
		}
		Q4TH(const std::vector<Node> &v)
		:Element<dim>(v)
		{
		    this->gaussor=2;
        }
		Q4TH(const std::vector<Node> &v,const unsigned int &go)
		:Element<dim>(v)
		{
		    this->gaussor=go;
        }
        inline ShapeFunctionGroup CreateShapeFunctionGroup()const
        {
        	ShapeFunctionGroup shfngr;

            //r,s,t,rs,st,rt,k

            shfngr.AddShapeFn(ShapeFunction( 0.25, 0.25, 0., 0.25, 0., 0., 0.25));//(1+r)(1+s)/4
            shfngr.AddShapeFn(ShapeFunction(-0.25, 0.25, 0.,-0.25, 0., 0., 0.25));//(1-r)(1+s)/4
            shfngr.AddShapeFn(ShapeFunction(-0.25,-0.25, 0., 0.25, 0., 0., 0.25));//(1-r)(1-s)/4
            shfngr.AddShapeFn(ShapeFunction( 0.25,-0.25, 0.,-0.25, 0., 0., 0.25));//(1+r)(1-s)/4

            return shfngr;
        }//Virtual
		inline const int DOF()const{ return 4; }	//TO MODIFY; THIS DEPENDS ON ELEMENT VARIABLES
		inline const int DIM()const{ return 2; }	//REDUNDANT
		inline const int DoFPerNode() const {return 1;}

		//If isotropic
		inline const Matrix<double> MatMatrix () const
		{
			cout << "Creating Element Matrix..."<<endl;
			Matrix<double> c(2,2);
			double k=this->mat->K(0);
			c[0][0] = c[1][1] = k;

			return c;


		}	//Create to save memory
		const unsigned int Field_Dim()const{return 1;}

	};


	//HEXA 8 Nodes
	template <int dim>
	class H8TH
	:public
	Element<dim>
	{

	public:
		H8TH():Element<dim>(){}
        H8TH(const std::vector<int> &v)
		:Element<dim>(v)
		{
            this->gaussor=2;
		}
		H8TH(const std::vector<int> &v, const unsigned int &go)
		:Element<dim>(v)
		{
            this->gaussor=go;
		}
		H8TH(const std::vector<Node> &v)
		:Element<dim>(v)
		{
		    this->gaussor=2;
        }
		H8TH(const std::vector<Node> &v,const unsigned int &go)
		:Element<dim>(v)
		{
		    this->gaussor=go;
        }
        inline ShapeFunctionGroup CreateShapeFunctionGroup()const
        {
            ShapeFunctionGroup shfngr;

            //Coeff Order
            //r,s,t,rs,st,rt,k, rst
            shfngr.AddShapeFn(ShapeFunction( 0.125, 0.125,-0.125, 0.125,-0.125,-0.125, 0.125,-0.125));//(1+r)(1+s)(1-t)/8
            shfngr.AddShapeFn(ShapeFunction(-0.125, 0.125,-0.125,-0.125,-0.125, 0.125, 0.125, 0.125));//(1-r)(1+s)/4
            shfngr.AddShapeFn(ShapeFunction(-0.125,-0.125,-0.125, 0.125, 0.125, 0.125, 0.125,-0.125));//(1-r)(1-s)/4
            shfngr.AddShapeFn(ShapeFunction( 0.125,-0.125,-0.125,-0.125, 0.125,-0.125, 0.125, 0.125));//(1+r)(1-s)/4
            shfngr.AddShapeFn(ShapeFunction( 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125));//(1+r)(1-s)/4
            shfngr.AddShapeFn(ShapeFunction(-0.125, 0.125, 0.125,-0.125, 0.125,-0.125, 0.125,-0.125));//(1+r)(1-s)/4
            shfngr.AddShapeFn(ShapeFunction(-0.125,-0.125, 0.125, 0.125,-0.125,-0.125, 0.125, 0.125));//(1+r)(1-s)/4
            shfngr.AddShapeFn(ShapeFunction( 0.125,-0.125, 0.125,-0.125,-0.125, 0.125, 0.125,-0.125));//(1+r)(1-s)/4

            return shfngr;

        }//Virtual
		inline const int DOF()const{ return 8; }	//TO MODIFY; THIS DEPENDS ON ELEMENT VARIABLES
		inline const int DIM()const{ return 3; }	//REDUNDANT
		inline const int DoFPerNode() const {return 1;}

		//If isotropic
		inline const Matrix<double> MatMatrix () const
		{
			cout << "Creating Element Matrix...HHH"<<endl;
			Matrix<double> c(3,3);
			double k=this->mat->K(0);
			c[0][0] = c[1][1] = c[2][2]=k;

			return c;


		}	//Create to save memory
		const unsigned int Field_Dim()const{return 1;}

	};

	template <int dim>
	//inline void Element < dim > ::Set_Jacobian()
	inline const GaussFullMatrices Element < dim > ::Jacobian()	const
	{
		GaussIntegrationScheme gi(this->GaussOrder(), dim);
		//Or elref.DIM()//
		GaussFullMatrices ret(dim, dim, gi);

		ShapeFunctionGroup shfngr = this->CreateShapeFunctionGroup();

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

		//this->shape_localgrad_matrices = temph;
		//this->jacobian = temph*this->XYZ();
		ret = temph*this->XYZ();

		return ret;
	}



	template <int dim>
	//TO MODIFY TO GaussMatrices
	inline const GaussFullMatrices Element<dim>::Shape_value_matrix() const
	{

		GaussFullMatrices m;

		GaussIntegrationScheme gi(this->GaussOrder(), dim);
		ShapeFunctionGroup shfngr = this->CreateShapeFunctionGroup();


		GaussFullMatrices ret(dim, shfngr.Size()*dim, gi);
		GaussFullMatrices retcomp(1, shfngr.Size(), gi);


		for (int g = 0; g < gi.NumPoints(); g++)
		{
			for (int h = 0; h < shfngr.Size(); h++)
			{
				//TO MODIFY
				retcomp[g][0][h] = shfngr.ShapeFn(h).Val(gi[g]);
				for (int comp = 0; comp < dim; comp++)
					ret[g][comp][(h * dim) + comp] = shfngr.ShapeFn(h).Val(gi[g]);
			}

		}

//		shape_value_matrices = ret;
//		shape_value_components=retcomp;
		return ret;
	}


} //FluxSol
#endif
