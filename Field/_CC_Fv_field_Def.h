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
#ifndef _CC_FVFIELD_DEF_H_
#define _CC_FVFIELD_DEF_H_

#include "Field.h"

namespace FluxSol
{

template <typename T>
_CC_Fv_Field<T>::_CC_Fv_Field(Fv_CC_Grid &grid)
{
	{		
		GridPtr=&grid;

		//Inicio variables PARA TODOS LOS CELLS (ver luego el borde
		for (int c=0;c<GridPtr->Num_Cells();c++)
		{
			T val;
			this->value.push_back(val);
		}

	_Fv_Boundary_Field <T> bf(grid.vBoundary());
	//Para igualar los boundaryfield debo ver el Scalar =
	BoundaryField=bf;

	}
}

//Must be verified the patch fields against mesh patches
template <typename T>
_CC_Fv_Field<T>::_CC_Fv_Field(Fv_CC_Grid &grid, const _Fv_Boundary_Field <T> &bfield)
{
	_Fv_Boundary_Field <T> bf = bfield;

	{		
		GridPtr=&grid;

		//Inicio variables PARA TODOS LOS CELLS (ver luego el borde
		for (int c=0;c<GridPtr->Num_Cells();c++)
		{
			T val;
			this->value.push_back(val);
		}

	//Para igualar los boundaryfield debo ver el Scalar =
	BoundaryField=bf;

	}
}
//
//template <typename T>
//void _CC_Fv_Field<T>::ToCellCenters(EqnSystem <T> &eqnsys)
//{
//	//Can compare number of cells vs number of eqn
//	for (int e=0;e<eqnsys.Num_Eqn();e++)
//	{
//		this->value[e]=eqnsys.Eqn(e).X();
//	
//	}
//
//}

//template<typename T>
//_CC_Fv_Field <T>::_CC_Fv_Field(InputFile &inputfile)
//{
//
//
//
//}


}

#endif

