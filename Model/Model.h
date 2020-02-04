#ifndef _MODEL_H
#define _MODEL_H

//#include "./FiniteVolume/Mesh/FvGrid.h"
//#include "./Solver/SIMPLE.h"
//#include "./Nastran/Nastran.h"
#include "./Nastran/Modelo.h"

using namespace FluxSol::Nastran;

namespace FluxSol{

//Modelo con Material
//Malla, archivo Nastran, BCs
//Pero sin analisis
class Model
{
	//En principio el archivo esta dentro del modelo
	//Archivo *nasfile;	//Cambiar el nombre de la clase a una que incluya el nombre Nastran
	Modelo *nasmodel;

protected:
	int bidim_model;


	//Devolucion de Miembros
	//Archivo & NasFile(){return (*nasfile);}

public:

	//Creo el 
	Model(){}


};

//class Fv_CC_Model:public Model
//{
//
//	Fv_CC_Grid *mesh;	//Recordar que el Campo tb linkea a la malla
//
//	public:
//	//Constructores
//	Fv_CC_Model(){};
//	Fv_CC_Model(const Fv_CC_Grid &);
//	Fv_CC_Model(const Archivo &nasfile, const Fv_CC_Grid &); //Genero un modelo a partir de un Nastran
//														 	 //Que a su vez genera la malla
//	
//	Fv_CC_Model(const Modelo &nasmod);
//	Fv_CC_Model(const std::string);
//	void Fv_CC_Model::Extract_Cells_and_BoundaryFromNastran();
//
//	////////////////////////////////////
//	//Funciones que devuelven miembros//
//	////////////////////////////////////
//
//	Fv_CC_Grid & Mesh(){return (*mesh);}
//
//
//};

};//Fin de FluxSol
#endif