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
#ifndef _VERTEX_H
#define _VERTEX_H


#include "../Type/Vec3d.h"
// ATENCION: ESTOS SON LOS EXTREMOS DE LAS CELDAS, ACA NO SE GUARDA INFORMACION
// DE PRESION NI VELOCIDAD
namespace FluxSol{
class _Vertex:public Vec3D{

    //vector <>
    private:
    vector<int> ind_cell_vecino;    //Cells que rodean al vertice
	
	//Parametros de nastran
	unsigned int SistCoord;				//Sistema de coordenadas
	unsigned int id;					//Id Nastran

    protected:

	public:
    _Vertex(double x=0., double y=0., double z=0.):Vec3D(x, y, z){}
	_Vertex(const int &i, const int &sc, vector<double> &coord):Vec3D(coord[0], coord[1], coord[2]){id=i;	SistCoord=sc;}




	};
}
#endif
