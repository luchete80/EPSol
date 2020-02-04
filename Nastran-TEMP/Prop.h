#ifndef PROP_H
#define PROP_H

#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <iomanip>
#include <locale>
#include <sstream>

using namespace std;

class Propiedad{
	public:
		Propiedad::Propiedad(){};
		Propiedad::Propiedad(const int &id, const int &mat);
		void Propiedad::Iniciar(const int &id, const int &mat);
		int Propiedad::VerId_Nastran();
		int Propiedad::Ver_Pos_Nastran();
	private:
		int Mat;
		int id;	//De nastran
		int pos_nastran;
};

#endif