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
#include "Cell.h"


using namespace std;

namespace FluxSol {

//Constructor con las conectividades
//_Cell::_Cell(vector <int> conect)
//{
//	_Cell();
//	id_vert=conect;
//}

//Inicio el cell con el id y las conectividades
//_Cell::_Cell(const int &Id, std::vector<int> conect)
//{
//	id=Id;
//	id_vert=conect;
//	this->num_vertex=conect.size();
//
//}


void _Cell::Conect(const vector<int> c)
{
    id_vert=c;
    num_vertex=c.size();
}

//
//string _Cell::Imprimir_Conect()
//{
//    ostream hola;
//    //Imprime una tira de enteros
//    for (int i=0;i<num_vertex;i++){
//
//    }
//}

string _Cell::Imprimir_Conect()//Devuelve un str de vertices
{

    string s,s2;
    char temp[100];
    for (int i=0;i<num_vertex;i++){
        sprintf(temp,"%d",Id_Vert(i));s2=temp;
        s+=s2+"  ";}
    return s;
}

//This is obtained from free-cfd solver
vector <int>  _Cell::GlobalVertFace(const int &nf)
{
	vector <int> tempNodes;
	int numfaceverts;
	//Recordar horario desde afuera segun numera NASTRAN
	int hexaFaces[6][4]= {
	{0,3,2,1},
	{4,5,6,7},
	{1,2,6,5},
	{0,4,7,3},
	{1,5,4,0},
	{2,3,7,6}
	};

	int prismFaces[5][4]= {
		{0,2,1,0},
		{3,4,5,0},
		{0,3,5,2},
		{1,2,5,4},
		{0,1,4,3}
	};
	int pyraFaces[5][4]= {
		{0,3,2,1},
		{0,1,4,0},
  		{1,2,4,0},
  		{3,4,2,0},
  		{0,4,3,0}
	};
	int tetraFaces[4][3]= {
		{0,2,1},
		{1,2,3},
  		{0,3,2},
  		{0,1,3}
	};

	//si ya se los switch no van
	switch (this->Num_Vertex()) {
		case 4: // Tetrahedra
			numfaceverts=3;
			tempNodes.assign(3,0);
			break;
		case 5: // Pyramid
			if (nf<1) {	//La primer cara tiene 4 lados (base) las restantes 3
				numfaceverts=4;
				tempNodes.assign(4,0);;
			} else {
				numfaceverts=3;
				tempNodes.assign(3,0);;
			}
			break;
		case 6: // Prism
			if (nf<2) { //La primer cara tiene 4 lados (base) las restantes 3
				numfaceverts=3;
				tempNodes.assign(3,0);
			} else {
				numfaceverts=4;
				tempNodes.assign(4,0);;
			}
			break;
		case 8: // Brick 
			numfaceverts=4;
			tempNodes.assign(4,0);;
			break;
	}//Fin del switch de cantidad de vertices por cell

	//esto lo puedo sacar
	//Asigno los vertices globales, que es tempNodes
	for (int fn=0;fn<numfaceverts;++fn) {
	switch (this->Num_Vertex()) {
		case 4: tempNodes[fn]=this->Id_Vert(tetraFaces[nf][fn]); break;
		case 5: tempNodes[fn]=this->Id_Vert(pyraFaces[nf][fn]);  break;
		case 6: tempNodes[fn]=this->Id_Vert(prismFaces[nf][fn]); break;
		case 8: tempNodes[fn]=this->Id_Vert(hexaFaces[nf][fn]);  break;
		} //Fin del switch
	}//End of for

	return tempNodes;

}

}