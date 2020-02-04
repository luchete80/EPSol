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
#ifndef _PATCHFIELD_H_
#define _PATCHFIELD_H_

#include "Field.h"

namespace FluxSol
{

template <typename T>
_Fv_Patch_Field<T>::_Fv_Patch_Field(const Patch &p)	
{
	Patch patch=p;

	this->numberofvals=0;

	for (int nf=0;nf<patch.Num_Faces();nf++)
	{
		T v(0.);
		this->value.push_back(v);
		this->numberofvals++;
	}
}


template <typename T>
void _Fv_Patch_Field<T>::AssignValue(const T &val)	
{
	for (int n=0;n<this->numberofvals;n++)
		this->value[n]=val;
}

}

#endif
