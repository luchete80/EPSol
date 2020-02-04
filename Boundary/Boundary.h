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
#ifndef _BOUNDARY_H
#define _BOUNDARY_H

#include <vector>
#include <list>
#include <string>
//#include "../Nastran/Varios.h"

namespace FluxSol{

//Varios patch forman un boundary
//CAN BE A TEMPLATE OR
//HERE face can refer either a point, line or face
//according to mesh dimension
class Patch{

    int id;
	int num_cells;
	int num_faces;

	std::vector <int> id_face;
	std::vector <int> id_cell;	//Este puede ddeducirse del anterior pero no la inversa

	const std::string name;


	public:
		Patch(){};
		Patch(const int &Id){id=Id;num_cells=num_faces=0;}
		Patch(const std::vector <int> cell, const std::vector <int> face){};
		//template <typename T,unsigned S> Patch( T (&v)[S]);
		Patch(std::list <int> &idfaces);

		void AddFace(const int &f){id_face.push_back(f);id_face[num_faces]=f;num_faces++;}
		int & Num_Faces(){return num_faces;};
		int & Id_Face(const int &i){return this->id_face[i];}


		// OJO QUE ESTA FUNCION DEBE CREARSE, VER EN DONDE
		// LO PIDE EL COMPILADOR
		Patch & operator=(const Patch&){return *this;}


	};

class Boundary{

    std::vector<Patch> patch;
	std::list <int> idpatch;	//Esta se recorre en orden
	std::vector <int> e;

	int num_faces;			//En la construccion se suma todo
	int num_patches;


    public:
    Boundary(){};
	//Boundary(std::vector<Patch> vp){patch=vp;}		Esto no lo puedo hacer porque no tengo un operador de patches
	Boundary(std::vector<int> vp){e=vp;}
	void AddPatch(Patch p){patch.push_back(p);}
	Boundary (std::vector <Patch> patchlist);
	int & Num_Faces(){return num_faces;}
	int & Num_Patches(){return num_patches;}
	Patch & vPatch(const int &i){return patch[i];}
};

}

#endif
