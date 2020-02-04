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

#ifndef _FIELDOPERATIONS_H_
#define _FIELDOPERATIONS_H_

#include "Field.h"
#include "../Type/Vec3D.h"
#include "../Type/Scalar.h"

namespace FluxSol
{


//THIS MUST BE DONE ONCE
//Inner product
template<typename T>
const _Field<T> & _Field<T>::operator &(const _Field<T> &right) const
{
	_Field<T> ret(this->Numberofvals());
	//Sizes must be equal and rank must be large than zero?
	for (int c = 0; c < this->Numberofvals; c++)
	{
		ret[c] = c&right[c];
	}

	return ret;
}


}


#endif
