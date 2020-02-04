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
#include "./Mesh/Grid.h"
#include "./Mesh/Structured.h"

using namespace FluxSol;

class _Rectangle_Grid:public _Grid{
private:
int num_cells_x,num_cells_y,num_cells_z;

public:
_Rectangle_Grid(){};
_Rectangle_Grid(const int &nex, const int &ney){num_cells_x=nex;num_cells_y=nex;};
};
