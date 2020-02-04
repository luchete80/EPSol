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
#ifndef _GRID_H
#define _GRID_H

#include "Cell.h"
#include "Face.h"
#include "./Mesh/Vertex.h"
#include "./Type/Vec3d.h"
#include "./Mesh/Node.h"
#include "Boundary.h"
#include <map>
#include <set>
#include <string>
#include <vector>

// ATENCION: ESTOS SON LOS EXTREMOS DE LAS CELDAS, ACA NO SE GUARDA INFORMACION
// DE PRESION NI VELOCIDAD

using namespace std;

namespace FluxSol{

//Like in free-cfd, a rawdata for mesh reading

class GridRawData { // This data will be destroyed after processing
public:
	std::vector<Vec3D> node;
	std::vector<int> cellConnIndex,cellConnectivity;
	std::vector< set<int> > bocoNodes; // Node list for each boundary condition region
	std::map<string,int> bocoNameMap;

	//Modificaciones del original free-cfd asumiendo en principio una sola zona
	//Modification of the original free-cfd
	std::vector<set<int> > bc_elem_list;	//Asummed at first one single zone
};

//Malla GENERICA
//de elementos finitos, o de volumenes finitos en cualquier variante
class _Grid{

    private:


    protected:

	GridRawData raw;
	//Archivos y logs
	//Que conviene???
	//std::ostringstream hola;
    //std::ofstream meshlog;


	bool escribir_log;				//Si escribir o no el log de creacion de malla
														//Podria hacerlo por niveles. En principio es todo o nada

    //Flags de iniciacion
    bool inicie_cells;
    bool inicie_faces;
	bool inicie_nodes;
	bool bidim_mesh;



	//Geometria
	//En principio esto puede ser elemento, o celda de cualquier tipo
    //Esta malla BASE solo tiene vertices y lista de nodos, no tiene nada mas
	std::vector <_Vertex> vert;			//Vertice, no es NI NODO NI ELEMENTO

	//LAS OTRAS COSAS DEPENDEN DEL TIPO DE MALLA

	int num_verts;						//Este Numero no es constante por si modifico la malla
														//Si voy a tener vertices tipo baffle no lo voy a tener en cuenta
	int num_cells_int;				//Numero de celdas efectivas
	int num_faces_int;				//Numero de faces efectivas (INTERIORES)
	int num_faces;
	int num_cells;

	int num_boundary_faces;

	//Para malla estructurada solamente
    vector <int> num_cells_xyz;
    vector <int> num_verts_xyz;

	string fileName;

	//Iteradores
	std::vector<_Vertex>::iterator vertit;
	std::vector<Node>::iterator nodeit;
	std::vector<_Cell>::iterator cellit;
    std::vector<_Face>::iterator faceit;

	//TO INHERITED MEMBERS
	std::vector <_Cell> cell;
	//This must be changed in FluxSol
	std::vector <Node>    node;    //En el centro de cada una de las cells



	//METODOS PRIVADOS -- NO CONVIENE HACERLOS PUBLICOS
	void Num_Verts(const int &nv);
	void Num_Verts_xyz(const int &x,const int &y,const int &z){
	num_verts_xyz.push_back(x);
	num_verts_xyz.push_back(y);
	num_verts_xyz.push_back(z);};
	void Num_Cells(const int &nc){num_cells=nc;};

	Boundary boundary;

	public:
	_Grid(){};
	//Assuming cgns
	_Grid(const char *name):fileName(name){/*this->Read_CGNS();*/};
	_Grid(const int &nex, const int &ney);
    _Grid(const int &nex, const int &ney, const int &lex, const int &ley);		//Uno de los tantos constructores, malla cuadrada de 4 nodos por elemento
    //formerly virtual
	void Iniciar_Caras();
    void Inicie_Cells(bool b){inicie_cells=true;}

	//Construct from cells
	_Grid(std::vector<_Cell>);

	//In those case with built in types, const return value is not importat
    int & Num_Verts(){return num_verts;}    //polimorfismo
	int & Num_Faces(){return num_faces;}    //polimorfismo
	int & Num_Boundary_Faces(){return num_boundary_faces;}    //polimorfismo
	int & Num_Cells() {return num_cells;}	//es el numero total
	void Bidim_Mesh (bool b){bidim_mesh=b;}
	std::vector <int> & Num_Verts_xyz(){return num_verts_xyz;}    //polimorfismo

    //Puedo hacer una funcion con varios argumentos en lugar de un vector para los indices
    vector <_Vertex> & Verts(vector<int> indvert);
	//Devuelvo el vertice v
	_Vertex & Vertex(const int &v){return vert[v];}



    //////////
    // LOGS //
    //////////

    //Manejo de vertices
    void PushVert(const _Vertex &v){vert.push_back(v);}
	void Insert_VertexVector( const std::vector<_Vertex> & vv){vert=vv;}

    //El log general no tiene mucho sentido que sea general (virtual)
	//formerly virtual
	void Log();
	void Trace();

    //void Log_Abrir(std::string nombre){meshlog.open(nombre.c_str());}
    void Log_Verts();

	void Read_InitialCGNS();
	//virtual void Read_CGNS();

	~_Grid(){}

};//_Grid

//Templates
template<typename T>
void Insertar_Cell(std::vector <int> con, T& t)
{
    t.cell.push_back(t);
}

}//Fin FluxSol

#endif
