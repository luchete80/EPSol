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
#include "Node.h"

// Constructor

FluxSol::Node::Node(const int &i, const int &sc, vector<double> &coord)
{
	id=i;
	Coord=coord;
	SistCoord=sc;
}

FluxSol::Node::Node(const int &i, const int &sc, const double &x, const double &y, const double &z)
{
	id = i;
	Coord.push_back(x);Coord.push_back(y);	Coord.push_back(z);
	SistCoord = sc;
}

//FluxSol::Node::Node(const int &i, const double &x, const double &y, const double &z)
//{
//	id = i;
//	Coord.push_back(x); Coord.push_back(y);	Coord.push_back(z);
//	SistCoord = 0;
//}


void FluxSol::Node::Node_Cargar_Ubic(const int u)
{
	pos_nastran=u;
}

