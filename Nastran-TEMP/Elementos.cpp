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
#include "Elementos.h"
#include "Varios.h"

//CBush::CBush()
//{}

Elemento::Elemento()
{
    es_cbush=false;

}

CBush::CBush()
{}

CBush::CBush(const int &i, vector <int> nodos)
{
	id=i;
	Nodo=nodos;
}
void Elemento::Iniciar_Prop(const int &i, const int &pos)
{
		idprop=i;
		pos_idprop_nastran=pos;
}

void Elemento::Iniciar_Prop_Modif(const int &v)
{
		id_nastran_prop_modificada=v;
}

const int Elemento::Ver_Prop_Modif()
{
	return id_nastran_prop_modificada;
}

int Elemento::Ver_Pos_Nastran()
{
	return pos_nastran;
}

void Elemento::Pos_Nastran(const int &p)
{
	pos_nastran=p;
}


void Elemento::Iniciar_Nodo_Master(const int &n)
{
	Elemento::nodo_master=n;
}

const int Elemento::Nodo_Master()
{
	return Elemento::nodo_master;
}

const int Elemento::Pos_prop()
{
	return pos_idprop_nastran;
}

void CBush::Iniciar_Caract(const int &i)
{
	Caract=i;
}


const int Elemento::VerId()
{
	return id;
}


Rigid::Rigid()
{}

Rigid::Rigid(const int &i, vector <int> nodos,bool es_spider)
{
	id=i;
	Nodo=nodos;
	Spider=es_spider;
}

bool Rigid::Es_Spider()
{
	if (Spider)
		return true;
	else
		return false;
}

//Devuelve un vector de conectividad del elemento
vector <int> Elemento::Conect()
{
	return Nodo;
}

void CBush::Iniciar_Es_Estructural(bool es)
{
	es_estructural=es;
}
const bool CBush::Iniciar_Es_Estructural()
{
	return es_estructural;
}

vector <int> Elemento::VerConect_Int()
{
	vector <int> nodos;
	for (int i=0;i<Nodos_int.size();i++)
		nodos.push_back(Nodos_int[i]);

	return nodos;
}

vector <int> Elemento::Conect_int()
{
	vector <int> v;
	for (int i=0;i<2;i++)
	v.push_back(Nodos_int[i]);
	return v;
}

void Elemento::Asignar_Conect_Interna(vector <int> v)
{
	Nodos_int=v;
}

void Elemento::Asignar_Prop(const int &id)
{
	pid=id;
}

const int Elemento::Linea_Nastran()
{
	return linea_nastran;
}
void Elemento::Iniciar_Fila_Nastran(const int &f)
{
	linea_nastran=f;
}

void Elemento::Asignar_Conect_Nastran(vector <int> v)
{
	Nodo=v;
}

void CBush::IniciarOri(int i)
{
	Ori=i;
}
const int CBush::VerOri()
{
	return Ori;
}

Elemento* Elemento::This()
{
	return this;
}

const int CBush::Ver_Caract()
{
	return Caract;
}

//Devuelve el indice INTERNO DE LA PROPIEDAD
const int Elemento::VerPid()
{
	return pid;
}

void Elemento::Iniciar_Nastran(const int &i, vector <int> nodos_nastran)	//Es con los nodos de nastran
{
	id=i;
	Nodo=nodos_nastran;
}

void Elemento::Leer_String(const std::string cad)
{
    //Primero miro todos los elementos que tienen los campos en el primero y en el segundo lugar (campos 4 y 5)
    if ( (cad.find("CBUSH",0)!=cad.npos) || (cad.find("CQUADR",0)!=cad.npos) || (cad.find("CQUAD4",0)!=cad.npos) || (cad.find("CBEAM",0)!=cad.npos) ||
         (cad.find("CBAR",0)!=cad.npos) || (cad.find("CTRIAR",0)!=cad.npos) || (cad.find("CTRIA3",0)!=cad.npos) )
    {
        //El primer campo siempre es el 4 en este caso
        Nodo.push_back(Leer_Campo(cad,4));//Es el indice de Nastran
        Nodo.push_back(Leer_Campo(cad,5));//Es el indice de Nastran
    }

    //RBE2 lee el campo 3 y 5
    if ( (cad.find("RBE2",0)!=cad.npos))
    {
        //El primer campo siempre es el 4 en este caso
        Nodo.push_back(Leer_Campo(cad,3));//Es el indice de Nastran
        Nodo.push_back(Leer_Campo(cad,5));//Es el indice de Nastran
    }

    //Si es CTRIA3leo un nodo mas
    if ( (cad.find("CTRIAR",0)!=cad.npos) || (cad.find("CTRIA3",0)!=cad.npos))
    {
        //El primer campo siempre es el 4 en este caso
        Nodo.push_back(Leer_Campo(cad,6));//Es el indice de Nastran
    }

    //Si es quad leo 2 nodos mas
    if ( (cad.find("CQUADR",0)!=cad.npos) || (cad.find("CQUAD4",0)!=cad.npos))
    {
        //El primer campo siempre es el 4 en este caso
        Nodo.push_back(Leer_Campo(cad,6));//Es el indice de Nastran
        Nodo.push_back(Leer_Campo(cad,7));//Es el indice de Nastran

    }

    if ( (cad.find("CBUSH",0)!=cad.npos))
    {
        Nodo.push_back(Leer_Campo(cad,3));//Es el indice de Nastran
    }

    id=Leer_Campo(cad,2);//Es el indice de Nastran

    if ( (cad.find("CBUSH",0)!=cad.npos))
        es_cbush=true;

	else  if ( (cad.find("RBE2",0)!=cad.npos))
    tipo="RBE2";

	//tipo=Leer_Hasta_Caracter(cad, ' ');
	//Leo propiedad
	pid=Leer_Campo(cad,3);


	//Coloco los primeros 8 caracteres de la cadena
	

}

void Elemento::Leer_StringVec(const std::vector <string> vcad)
{

    bool es_rigid=false;

    for (int i=0;i<vcad.size();i++)
    {
        string cad=vcad[i];
        if (i==0)
        {
            //Primero miro todos los elementos que tienen los campos en el primero y en el segundo lugar (campos 4 y 5)
            if ( (cad.find("CBUSH",0)!=cad.npos) || (cad.find("CQUADR",0)!=cad.npos) || (cad.find("CQUAD4",0)!=cad.npos) || (cad.find("CBEAM",0)!=cad.npos) ||
                 (cad.find("CBAR",0)!=cad.npos) || (cad.find("CTRIAR",0)!=cad.npos) || (cad.find("CTRIA3",0)!=cad.npos) )
            {
                //El primer campo siempre es el 4 en este caso
                Nodo.push_back(Leer_Campo(cad,4));//Es el indice de Nastran
                Nodo.push_back(Leer_Campo(cad,5));//Es el indice de Nastran
            }

            //RBE2 lee el campo 3 y 5
            if ( (cad.find("RBE2",0)!=cad.npos))
            {
                //El primer campo siempre es el 4 en este caso
                Nodo.push_back(Leer_Campo(cad,3));//Es el indice de Nastran
                Nodo.push_back(Leer_Campo(cad,5));//Es el indice de Nastran

                int campo=6;
                int val;
                while(Leer_Campo(cad,campo)!=0)
                {
                    val=Leer_Campo(cad,campo);
                    Nodo.push_back(val);
                    campo++;
                }
            }

            //Si es CTRIA3leo un nodo mas
            if ( (cad.find("CTRIAR",0)!=cad.npos) || (cad.find("CTRIA3",0)!=cad.npos))
            {
                //El primer campo siempre es el 4 en este caso
                Nodo.push_back(Leer_Campo(cad,6));//Es el indice de Nastran
            }

            //Si es quad leo 2 nodos mas
            if ( (cad.find("CQUADR",0)!=cad.npos) || (cad.find("CQUAD4",0)!=cad.npos))
            {
                //El primer campo siempre es el 4 en este caso
                Nodo.push_back(Leer_Campo(cad,6));//Es el indice de Nastran
                Nodo.push_back(Leer_Campo(cad,7));//Es el indice de Nastran

            }

            id=Leer_Campo(cad,2);//Es el indice de Nastran

            if ( (cad.find("CBUSH",0)!=cad.npos))
                es_cbush=true;

            if ( (cad.find("RBE2",0)!=cad.npos))
            {
                tipo="RBE2";
                es_rigid=true;
            }


        }
        else    //No es la primera linea
        {
            if (es_rigid)   //Leo los slaves de mas
            {
                int campo=2;
                int val;
                while(Leer_Campo(cad,campo)!=0)
                {
                    val=Leer_Campo(cad,campo);
                    Nodo.push_back(val);
                    campo++;
                }

            }


        }

    }
	

}


