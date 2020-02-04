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
#ifndef _FIELDDEF_H_
#define _FIELDDEF_H_

#include "Field.h"


using namespace std;
namespace FluxSol {


//template <typename T>
//_Field<T>::Log ()
//{
	//for ()

//}
//
//template <typename T>
//Vertex_Fv_Field<T>::Vertex_Fv_Field(Fv_CC_Grid &grid)
//{
//	GridPtr=&grid;
//
//	T val;
//	int nv= grid.Num_Verts();
//	value.assign(nv,val);	
//}

template <typename T>
void _Field<T>::operator=(const T &val)
{
	for (int v=0;v<value.size();v++)
	{
		value[v]=val;
	}
}

//template <typename T>
//void _Field<T>::operator=(const double &val)
//{
//	for (int v=0;v<value.size();v++)
//	{
//		value[v]=val;
//	}
//}

template <typename T>
_Field<T> _Field<T>::operator=(const double &val)
{
	for (int v=0;v<value.size();v++)
	{
		value[v]=val;
	}
	return *this;
}

template <typename T>
_Field<T>::_Field(const int &numval, const double &number =0.)
{
	for (int v=0;v<numval;v++)
	{
		T num(number);
		value.push_back(num);
	}
}


template <typename T>
const Scalar MaxDiff(const _Field<T> &field1, const _Field<T> &field2)
{
	Scalar r;
	//Check for 
	if (field1.Numberofvals() != field2.Numberofvals() )
		return r; 

	Scalar maxdiff(1.0e10);
	for (int c=0;c<field1.Numberofvals();c++)
	{
		//Can be used the general Norm function
		Scalar diff=(field1[c].Norm()-field2[c].Norm())/field1[c].Norm();
		//Can define an operator
		if (diff<maxdiff)
			maxdiff=diff.Norm();
	}
}



// INITIALLY NO VIRTUAL FUNCTION
//template <typename T>
//_Field<T> & _Field<T>::operator&(const _Field<T> &right) const
//{


//}

//Inner Product Field
template<typename T>
SurfaceField<typename innerProduct < T, T> ::type> operator &(const SurfaceField<T> &left,const SurfaceField<T> &right)
{
	SurfaceField<typename innerProduct < T, T> ::type> ret(left.Numberofvals());
	typename innerProduct < T, T> ::type val;
	//Sizes must be equal and rank must be large than zero?
	for (int c = 0; c < left.Numberofvals(); c++)
	{
		val = left.Val(c) & right.Val(c);
		ret.Val(c,val);
	}

	return ret;
}



}//FluxSol


#endif
