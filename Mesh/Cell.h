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
#ifndef _CELL_H_
#define _CELL_H_

#include <vector>
#include <string>
#include <stdio.h>
#include "../Mesh/Node.h"
#include "../Type/Scalar.h"

namespace FluxSol {
class _Cell
{
    private:

    protected:


	//const string type;							//Hexa, Tetra,Prism,
	string type;							//Hexa, Tetra,Prism,
	std::vector <int> id_neigbourcell;			// Celdas vecinas
	std::vector <int> id_face;						// Cara
	std::vector <int> id_vert;                       //Vertices NO SON NODOS
    //El indice del nodo podria ser el mismo que el de cell no/
    std::vector <int> id_node;                       //Ojo aca cambia
													 //Estos son stencils??
	Scalar vp;							//Cell volume
	int num_faces;
	int num_vertex;
	int id;                             //El propio id

	bool es_coordenada; //Para evitar correccion por no ortogonalidad

	public:
		//Constructores
		_Cell(){};
		_Cell(std::vector <int> conect)
		{	
		this->id_vert=conect;
		}				//Vector de conectividades en forma de lista
		_Cell(const int &Id, std::vector<int> conect)
		{
			this->id=Id;
			id_vert=conect;
			this->num_vertex=conect.size();
		}

		// VARIABLES DE CEDAS
		int &  Num_Vertex(){return num_vertex;}
		int &  Num_Faces (){return num_faces; }

		void Conect(const std::vector<int> c);
		std::vector<int> & Conect(){return id_vert;}
		int & Vert(const int &i){return id_vert[i];}

        int Id(){return id;};   //Devuelve el id del Cell
		void Id(const int &Id){id=Id;}

		//Funcion para ver el id de un vertice
		int & Id_Vert (int i){return id_vert[i];}
		int & Id_Face (int i){return id_face[i];}
		void Id_Face (const int &i, const int &val){id_face[i]=val;}

		void Init_Idface(const int &i){id_face.assign(i,0);}
		void Init_NumberofFaces(const int &i){num_faces=i;}
		void Id_Vert(const int &i, int &pos){id_vert[pos]=i;}   //Conectividades por indice
        void Id_Vert(std::vector <int> i){id_vert=i;}           //Conectividades por vector

        //Log, cada una de las celdas tiene uno distinto
		void Log(){};
		std::string Imprimir_Conect();//Devuelve un str de vertices

		Scalar& Vp(){return vp;}

        //Para ver las variables de cell
        int Id_Node(int i){return id_node[i];}
        std::vector <int> Id_Node(){return id_node;}

		//const 


		//Faces
		void AddIdFace(const int &id){id_face.push_back(id);num_faces++;};
		vector <int> GlobalVertFace(const int &nf);
		
		//Neighbours
		void AssignNeighbours(std::vector <int> nc) {id_neigbourcell=nc;};
		void AddNeighbour(const int &i) {id_neigbourcell.push_back(i);};
		int Neigbour(const int &i){return id_neigbourcell[i];}
		unsigned int Num_Neighbours(){return id_neigbourcell.size();}

		~_Cell(){};


};

}//FluxSol

#endif
