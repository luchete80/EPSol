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
#include "Analysis.h"
using namespace FluxSol;
using namespace FluxSol::Nastran;

//Analysis::Analysis ( string archivo, const Fv_CC_Grid &mesh)
//{
//	// Abro el archivo 
//	//nasfile = new Archivo;
//	//nasfile->Iniciar(archivo);
//
//	//Ojo que los elementos de beam no son nuevos elementos (cells), sino 
//	//condiciones de borde
//	//Modelo para copiar las caracteristicas
//	Modelo modnas(archivo);	//Modelo de nastran
//
//	Fv_CC_Model fvmodel(archivo);	//Cargo el Fv model y dps le asigno su mesh a mesh
//	//mesh=fvmodel->mesh;
//	
//	//modnas.Convert2FvGrid();
//	//int numverts=nasfile->Buscar_numnodos();	//Ojo aca que los nodos no son asi
////	vector <Nodo> vert;
////	vert=nasfile->Leer_Nodos();
//
//}

