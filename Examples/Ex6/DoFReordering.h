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
#include "EPSol.h"


#ifndef _DOF_REORDERING_
#define _DOF_REORDERING_

using namespace FluxSol;
//INHERITING A CLASS TO TEST AND IMPLEMENT IT
template <int dim>
class DoFHandlerAndReordering:
public DoFHandler < dim >
{
	protected:
	
	public:

};


#endif
