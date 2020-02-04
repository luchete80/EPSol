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
#include "SistCoord.h"

SistCoord::SistCoord(int i, int tipo)
{
	id=i;
	Tipo=tipo;
}


const int SistCoord::VerTipo()
{
	return Tipo;
}

//Devuelve id de NASTRAN
const int SistCoord::VerId()
{
	return id;
}

//
//
//
//
//
//////Devuelve el tipo de sistema de coordenadas en funcion 
////Elvector Sc se divide en id y tipo
//




