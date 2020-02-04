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
#include "Material.h"


Material::~Material()
{}

const int Material::VerId()
{
	return id;
}

Material::Material(Caract c, const int &i,const double &h, const double &p, const double &area)
{
	E=1.0E08;
	nu=0.3;
	id=i;
	huelgo=h;
	//Defino con la caracteristica los puntos y el Et
	sy=c.VerMu()*c.VerPo()/area;	//s=muxpregarga
	Et=p*c.VerMu()*c.VerPo()/(huelgo);

}

const double Material::Ver_Sy()
{
	return sy;
}

const double Material::Ver_Et()
{
	return Et;
}

Material* Material::This()
{
	return this;
}

void Material::Crear_strings()
{
	//Creo string de material bilineal
	string cad="MAT1    ";
	cad+=intastring_alineado (id);
	cad+=dobleastring_alineado(E);
	cad+="      0.";	//G
	cad+=dobleastring_alineado(nu);
	cad+='\n';
	cad+="MATS1   ";
	cad+=intastring_alineado (id);
	cad+="        ";	//Tables1, no refiere a ninguna tabla
	cad+=" PLASTIC";
	double H=Et/(1.0-Et/E);
	cad+=dobleastring_alineado(H);
	cad+="       1";	//Fluencia por Von Mises
	cad+="       1";	//Endurecimiento isotropico
	cad+=dobleastring_alineado(sy);	//Endurecimiento isotropico
	cadbilineal=cad;
}

const string Material::Cadena_Bilineal()
{
	return cadbilineal;
}