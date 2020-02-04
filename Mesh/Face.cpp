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
#include "Face.h"

using namespace std;

namespace FluxSol{

_Face::_Face(const vector<int> & cells, const vector<int> & verts):null_flux_face(false)
{
	_Face();
    cell=cells;
    //Inicio los punteros de vertices
    //Esto es mucho mas facil para crear aca las superficies y otras cosas

}

_Face::_Face(const int &Id, const vector <_Vertex> &verts):null_flux_face(false)
{
    id=Id;   //Indice de la cara
    //Calculo posicion del baricentro
    sm=Prom(verts);

    af=CalcNormal(verts);
    norm_af=fabs(af);

	num_vertex=verts.size();

}


Vec3D &CalcNormal(const vector <_Vertex> &verts)
{
    Vec3D res;
	Vec3D ve;
    for (int v=3;v<verts.size();v++)
    {
		ve = verts[v - 1] - verts[0];
		//ve = verts[v - 1] - verts[0];
        //This works
		//res+= (ve).cross((verts[v]-verts[0]));
		//OpenFoam cross product sintax
		res += (verts[v - 1] - verts[0]) ^ (verts[v] - verts[0]);
		//res += (ve).cross((verts[v] - verts[0]));
		//res += (verts[v - 1] - verts[0]).cross((verts[v] - verts[0])); //PARAMETERS OF THIS FUNCTION MUST BE CONST 
																		//BINARY OPERATOR
		//res+= ve.cross((verts[v]-verts[0]));

	}

    res*=0.5;
    return res;
}

}
