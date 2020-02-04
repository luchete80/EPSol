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
#ifndef _FVFIELD_H
#define _FVFIELD_H

#include "Field.h"
#include "GeometricField.h"
//A general FvField can inherit fields below

namespace FluxSol
{

	template <typename T>
	class Vertex_Fv_Field :public _Field<T>
	{
		Fv_CC_Grid *GridPtr;

	public:
		Vertex_Fv_Field(){};
		Vertex_Fv_Field(Fv_CC_Grid &grid);

		//Later it will be inherited
		Fv_CC_Grid  & Grid(){ return *GridPtr; }

	};

	/////////////////////////////
	///// CAMPOS GEOMETRICO /////
	/////////////////////////////
	// IMPORTANTE REFERENCIA A LA MALLA
	//Tienen asociadas mallas, siempre se llama igual *GridPtr
	//y puede ser distinto tipo de malla
	//CONVIENE PONER FvField como padre

	//Campo de superficie
	template<typename T>
	class _Surf_Fv_Field :public _Field<T>{

		Fv_CC_Grid *GridPtr;       //Se podria probar con un puntero general
		// en la clase base


	protected:

	public:
		_Surf_Fv_Field(const int &numval, const double &value = 0.) :_Field(numval, value){}
		_Surf_Fv_Field(){};         //Constructor
		Fv_CC_Grid  & Grid(){ return *GridPtr; }
		//Adding boundary face
		_Surf_Fv_Field(Fv_CC_Grid &grid)
		{
			GridPtr = &grid;
			//LO HAGO PARA TODOS LOS FACES
			for (int f = 0; f<GridPtr->Num_Faces(); f++)
			{
				if (GridPtr->Face(f).Boundaryface() && !GridPtr->Face(f).Is_Null_Flux_Face())
				{
					T val;
					this->value.push_back(val);
					this->idgrid_domain.push_back(f);
				}
			}

		}

		//Devuelve el valor en una Surface 
		T& GridIdSurFaceVal(const int &i)		//El valor de una 
		{
			//Tengo que ver cual es la face i del field
			for (int fid = 0; fid<this->Numberofvals(); fid++)
			{
				if (i == this->Idgrid_domain(fid))
					return this->value[fid];
			}
		}

	};

	//boundary type

	enum PatchFieldType { FIXEDVALUE, FIXEDGRADIENT };

	//
	template<typename T>
	class _Fv_Patch_Field :
		public _Surf_Fv_Field<T>
	{
	protected:
		int patchid;
		PatchFieldType type;


	public:
		_Fv_Patch_Field(){};
		explicit _Fv_Patch_Field(const Patch &);
		int & PatchId(){ return patchid; };
		//This must be in parent class
		void AssignValue(const T &val);

		PatchFieldType& Type(){ return type; }


	};

	//Diriclet Field
	template <typename T>
	class _Fv_FixedValue_Patch_Field :public _Fv_Patch_Field<T>
	{


	public:
		explicit _Fv_FixedValue_Patch_Field(const Patch &p) :_Fv_Patch_Field(p){ type = FIXEDVALUE; };
	};

	//Newmann Field
	template <typename T>
	class _Fv_FixedGradient_Patch_Field :public _Fv_Patch_Field<T>
	{


	public:
		explicit _Fv_FixedGradient_Patch_Field(const Patch &p) :_Fv_Patch_Field(p){ type = FIXEDGRADIENT; };
	};


	//Instead a FieldField as in OpenFoam has a vector in Fields
	//This will be soon a field
	template<typename T>
	class _Fv_Boundary_Field
	{
	protected:
		//Can be Fixed Value and Fixed Gradient
		std::vector < _Fv_Patch_Field < T > > patchfield;

	public:
		//_Fv_Boundary_Field(const _CC_Fv_Field < T > &);
		_Fv_Boundary_Field(){}
		_Fv_Boundary_Field(const Boundary &bound)
		{
			Boundary b = bound;
			//Generate a patchfield for every patch in boundary
			for (int np = 0; np<b.Num_Patches(); np++)
			{
				_Fv_Patch_Field < T > pf(b.vPatch(np));
				patchfield.push_back(pf);
			}

		}

		_Fv_Patch_Field < T > & PatchField(const int &i){ return this->patchfield[i]; }
		void AddPatchField(_Fv_Patch_Field<T> &field)  { patchfield.push_back(field); }

	};

	//CAMPO VOLUMETRICO DE VOLUMENES FINITOS
	//CONVIENE QUE DERIVE 
	// La T es Scalar, vector o tensor
	template<typename T>
	class _CC_Fv_Field :public _Field<T>{

		//Busco en los nodos de GridPtr
		//El vector value refiere a cada uno de los node

		//

		//SET DE DIMENSIONES

	protected:
		//Conviene colocar referencia o puntero??
		_Fv_Boundary_Field <T> BoundaryField;	//En las caras


	public:

		//EVALUATE IF IT IS CONST
		Fv_CC_Grid  *GridPtr;       //Se podria probar con un puntero general	

		_CC_Fv_Field(Fv_CC_Grid &grid);
		_CC_Fv_Field(Fv_CC_Grid &grid, const _Fv_Boundary_Field <T> &bfield);

		Fv_CC_Grid  & Grid(){ return *GridPtr; }

		const Fv_CC_Grid  & ConstGrid()const { return *GridPtr; }

		//_CC_Fv_Field (InputFile &inputfile);

		_CC_Fv_Field(){};

		//Conveccion, interpolacion o reconstruccion, de solo una celda
		//THIS IS OLD. FIELD DOES NOT RETURN AN EQ SYSTEM BECAUSE REFER TO EQSYS AND THIS TO FVGRID
		//AT THE SAME TIME FVGRID REFERS TO GEOMFIELD AND THIS TO FIELD, AND FIELD TO EQSYS
		//const Eqn <T> ToFaces(const int &icell, const std::vector <Scalar> flux) const;
		
		_Fv_Boundary_Field <T> &  Boundaryfield(){ return BoundaryField; }

		//void ToCellCenters(EqnSystem <T> &eqnsys);


		_CC_Fv_Field<T> operator=(const double &val)
		{
			for (int v = 0; v<value.size(); v++)
			{
				value[v] = val;
			}
			return *this;
		}


	};


	////////////////////////////
	// FUNCIONES DE TEMPLATES //
	////////////////////////////
	//Reconstruccion a las caras
	//Esta funcion es llamada por la divergencia
	//Devuelve una ecuacion que es como un stencil pero tiene a la vez la componente
	//de source
	//Esta funcion devuelve Sum n (-min[mf,0])
	//O lo que es lo mismo Si m<0 fi f=fi p si no fif=fiN
	//template < typename T >
	//const Eqn <T> _CC_Fv_Field<T>::ToFaces(const int &icell, std::vector <Scalar> flux) const
	//{
	//	vector <T> an;      //coeficientes de los vecinos
	//	T ap(0.);          //Lo inicio vacio, podria tb igualarlo a 0
	//	//Recorro las caras del cell
	//	//En el caso de upwind o CDS simple la molecula es siempre la misma

	//	for (int iface = 0; iface<GridPtr->Cell(icell).Num_Faces(); iface++)
	//	{
	//		T coef(0.);
	//		an.push_back(coef); //Agrego un vecino mas
	//		//
	//		if (flux[iface].Val()>0.)   //Flujo saliente, va a ap
	//			ap += flux[iface].Val();       //Este oprador suma un valor a derecha constante

	//		else if (flux[iface].Val()<0.)  //Flujo entrante, va al vecino
	//			an[iface] += flux[iface].Val();

	//	}
	//	//Creo la ecuacion con el coeficiente central y los vecinos
	//	//Como no especifico la fuente, esta es nula
	//	Eqn <T> eqn(ap, an);
	//	return eqn;


	//	//Field Output
	//	void Write_vtk();
	//}

	template <typename T>
	const Scalar MaxDiff(const _Field<T> &field1, const _Field<T> &field2);
	//#include "Field_Def.h"

}//FluxSol

#endif
