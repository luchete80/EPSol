#ifndef NODO_H
#define NODO_H

#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>

using namespace std;

class Nodo
{
public:
	Nodo();
	Nodo(const int &i, const int &sc, vector<double> &coord);
	Nodo(const double &d);
	int VerId_Nastran();
	vector <double> Coords();
	const int Sc();
	const int Sc_int();
	const int Pos_Nastran();
	void Nodo_Cargar_Ubic(const int u);
	void Cargar_IdSc_Int(const int &val);
	void Iniciar_Coord_cart(vector<double> &coord);
	vector <double> Coord_carts();

	private:
	int id;	//Id de Nastran
	std::vector <double> Coord;
	std::vector <double> Coord_cart;
	int SistCoord;
	int SistCoord_int;
	int pos_nastran;	//Ubicacion en el NASTRAN

};

#endif
