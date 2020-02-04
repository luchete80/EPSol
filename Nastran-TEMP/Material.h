#ifndef MATERIAL_H
#define MATERIAL_H

#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>

#include "./Nastran/Nodo.h"
#include "./Nastran/Elementos.h"
#include "./Nastran/Caract.h"
#include "./Nastran/Cadenas.h"

using namespace std;

class Material{

	public:
		Material::Material(){};
		Material::Material(Caract c, const int &i,const double &h, const double &p, const double &area);
		Material::~Material();
		const int VerId();
		Material* Material::This();
		void Material::Crear_strings();
		const string Material::Cadena_Bilineal();
        const double Material::Ver_Sy();
		const double Material::Ver_Et();

	private:
		string cadbilineal,cadmultilineal;
		double nu;		//Se supone elastoplástico con una 2da rigidez
		int id;	//Id en NASTRAN
		double E,Et;	//Et para bilineal
		double sy;
		double huelgo;
		string nastranbi,nastranmulti;

};


#endif
