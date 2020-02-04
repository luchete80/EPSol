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
#include "Nodo.h"

Nodo::Nodo()
{
}

Nodo::Nodo(const int &i, const int &sc, vector<double> &coord)
{
	id=i;
	Coord=coord;
	SistCoord=sc;
}

int Nodo::VerId_Nastran()
{
	return id;
}

const int Nodo::Sc()
{
	return SistCoord;
}

void Nodo::Iniciar_Coord_cart(vector<double> &coord)
{
	Coord_cart=coord;
}

const int Nodo::Sc_int()
{
	return SistCoord_int;
}


vector <double> Nodo::Coords()
{
	return Coord;
}

vector <double> Nodo::Coord_carts()
{
	return Coord_cart;
}

const int Nodo::Pos_Nastran()
{
	return pos_nastran;
}

void Nodo::Nodo_Cargar_Ubic(const int u)
{
	pos_nastran=u;
}

void Nodo::Cargar_IdSc_Int(const int &val)
{
	SistCoord_int=val;
}

