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
#include "Varios.h"

using namespace std;

//El campo empieza desde 1 !!!!!!!!!!!!!!!!!!!!
int Leer_Campo(const string cad, const int campo)
{
    string strlin=cad;
    string cadtemp;
    int iniciocad=(campo-1)*8;

    for (int k=0;k<8;k++)
    {
        cadtemp+=strlin[iniciocad+k];

    } //Delf for k

        int Number;
        //if ( ! (istringstream(cadtemp) >> Number) ) Number = 0;
        //Si no es un entero
        if ( (cadtemp.find(".",0)!=cadtemp.npos))           Number=0;

        return Number;
}

template <typename T>
void listar_vec(const std::vector <T> v, ofstream f)
{
    for (int i=0;i<v.size();i++)
    f<<v[i]<<endl;
}

int & SearchVal(const int &i,std::vector <int> &v)
{
	int r=-1;
	for (int vi=0;vi<v.size();vi++)
		if (v[vi]==i)
			r=vi;

	return r;
}

//Devuelve verdadero si los encontro a todos los ind en v
inline bool FindAllVals(std::vector <int> &ind, std::vector <int> &v)
{
	bool r=true;

	bool enc=true;
	int vi=0;
	while (vi<v.size() && enc)
	{
		bool enc_elem=false;
		for (int vind=0;vind<ind.size();vind++)
			if (ind[vind]==v[vi])
				enc=true;
		if (!enc_elem)
			enc=false;
		vi++;
	}
	return r;
}

//Indica la cantidad de elementos del vector ind que encuentra en v
//Si no encuentra devuelve 0
int & NumVecElemFound(std::vector <int> &ind, std::vector <int> &v)
{
	int r=0;

	int vi=0;
	while (vi<v.size())
	{
		bool enc_elem=false;
		for (int vind=0;vind<ind.size();vind++)
			if (ind[vind]==v[vi])
				r++;
		vi++;
	}
	return r;

}

//Valores en comun que tienen
std::vector <int> & CommonValues(std::vector <int> &ind, std::vector <int> &v)
{
	vector <int> r;

	int vi=0;
	while (vi<v.size())
	{
		bool enc_elem=false;
		for (int vind=0;vind<ind.size();vind++)
			if (ind[vind]==v[vi])
				r.push_back(v[vi]);
		vi++;
	}
	return r;
}

bool ValueFound(const int &val, std::vector<int> vec) 					
{	
	bool esta=false;
	for (int nvint=0;nvint<vec.size();nvint++)
	{
		if (vec[nvint]==val)
			esta=true;
	}
	return esta;

}
