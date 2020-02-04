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
#ifndef _SHAPE_H
#define _SHAPE_H

#include <vector>
#include "../Type/Vec3d.h"
#include "Function.h"

using namespace std;

namespace FluxSol
{
//Generic Geometric Shape
//Templatized?
class Shape{

protected:
    unsigned int dim;

    std::vector <int> vertex;                    //indice del vertice
	std::vector <int> cell;
	int num_vertex;
    //Directamente puedo colocar el puntero a cada uno de los vertices
	//Es un vector de punteros, mucho mas pequeño en memoria que los indices
	std::vector <Vec3D> *pvertex;
    Function <double> ds;       //Size diferential, function

public:
    //Constructors
    Shape(){}
    virtual void Set_ds(){};
    //Allows to construct from Vertex and Vertex Inherited classes (such as nodes)
    Shape(const int &Id,const std::vector <Vec3D *> &verts);

    const int & Cell(const int &i)const {return cell[i];};
	const int & Vert(const int &i)const {return vertex[i];};
    const int & NumCells(){return cell.size();}
    const int & NumVerts(){return num_vertex;}
    virtual const Function < double > & dS()const{};       //May be called as template too

};

}//FluxSol
#endif
