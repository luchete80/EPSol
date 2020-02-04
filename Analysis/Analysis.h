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
#ifndef _ANALYSIS_H_
#define _ANALYSIS_H_

#include "./Nastran/Nastran.h"
#include "./Nastran/Modelo.h"

#include "./Model/Model.h"

namespace FluxSol {

//En un principio entiendo que el analisis se genera a partir de Nastran
class Analysis
{
	//TRANQUILAMENTE PUEDE TENER OTROS ARCHIVOS QUE NO SEAN NASTRAN
	//Archivo *nasfile;	En principio el archivo esta dentro del modelo
	//Fv_CC_Grid *mesh;	//Recordar que el Campo tb linkea a la malla
                        //Pueden ser varios analysis heredados, o varios punteros

	//Solver


public:

	Analysis(){}
	//En este no me importa el modelo Nastran
	//Analysis(string, const Fv_CC_Grid &);


};

};

#endif
