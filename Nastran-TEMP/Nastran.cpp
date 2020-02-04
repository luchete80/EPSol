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
#include "Nastran.h"

using namespace FluxSol;

//Este constructor por defecto se llama cuando se crea desde otra clase!!!
Archivo::Archivo()
{
}

const vector<int> Archivo::Buscar_ini_fin(string cadena)
{
	vector <int> inifin;

	bool fin=false;
	int veces=0;
	int j=0;
	while(!fin && j<numfilas)
	{
	   //Veo si es una linea de definicion de nodos
	   //if(strncmp(lineas[j].c_str(),cadena.c_str(),4))
	     if(lineas[j].find(cadena,0)!=lineas[j].npos)
		   fin=true;
	   j=j+1;
	}

	if (fin)
		inifin.push_back(j-1);
	else
	{
		inifin.push_back(-1);
		inifin.push_back(-1);
		return inifin;	//No encontro nada
	}
	
	//Busco desde el final
	fin =false;
	j=lineas.size()-1;
	while(!fin && j>0)
	{
	   //Veo si es una linea de definicion de nodos
	   //if(strncmp(lineas[j].c_str(),cadena.c_str(),4))
		if(lineas[j].find(cadena,0)!=lineas[j].npos)
		   fin=true;
	   else
		   j=j-1;

	}

	if (fin)
		inifin.push_back(j);
	else
		inifin.push_back(-1);

	return inifin;

}


//Busca cuantas veces se repitio un comando o string
int Archivo::Buscar_num_com(string cadena)
{
	int veces=0;
	for (int i=1;i<numfilas;i++)
	{   //Veo si es una linea de definicion de nodos
	   //if(strncmp(lineas[j],cadena,4))
		if(lineas[i].find(cadena,0))
		   veces=veces+1;
	}
	return veces;
}

//Constructor que abre el archivo
void Archivo::Iniciar(string cad)
{
	nombre=cad;
	ifstream fen (cad.c_str(), ifstream::in);
	cout << "Creando Cadenas (Si demora modificar push_back).....) " <<endl;
	//Abro archivo
	cout<< "Abriendo archivo..."<<cad<<endl ;
	int i=0;
	string s;
	while(getline(fen, s)) // Discards newline char
	{
		lineas.push_back(s);
		i=i+1;
	}

	//Escribo el vector posicion de nodos
	if (i!=0)
	{
		numfilas=i;
		Buscar_pos_nodos();
	}
	else	//archivo vacío
	{
		cout << "Nastran file is empty or not found ..." <<endl;
	}
}

void Archivo::Buscar_pos_nodos()
{
	pos_nodos=Buscar_ini_fin("GRID    ");
}


//Devuelvo un vector de nodos
//Ya debe haber iniciado
const vector <Nodo> Archivo::Leer_Nodos()
{
	vector <Nodo> Nodos;
	Nodo nini;
	int id,idsc;
	//REVISAR ESTO!!!
	int iniciocad[]={24,32,40};
	vector <bool> modifcoord;
	vector <double> vectemp;
	vectemp.assign(3,0.0);
	modifcoord.assign(3,false);
	int nnodos=pos_nodos[1]-pos_nodos[0]+1;
	Nodos.assign(nnodos,nini);

	//ASUMO QUE LOS NODOS SIEMPRE TIENEN LAS COSAS ALINEADAS
	cout<<"Leyendo Nodos"<<endl;
	int indn=0;
	for (int i=pos_nodos[0];i<pos_nodos[1]+1;i++)
	{
		sscanf(lineas[i].c_str(),"%*s %d %d",&id,&idsc);
		string strlin=lineas[i];

		for (int dim=0;dim<3;dim++)
		{
			string cadtemp="";
			//Leo de a uno los 8 caracteres de formato NASTRAN
			int inc=0; //La cadena que leo puede tener mas caracteres
					//que 8 si hay que agregarle el E
			for (int k=0;k<8;k++)
			{
				cadtemp+=strlin[iniciocad[dim]+k];
				//Chequeo si el - está dps del punto
				//if (cadtemp.find("-",0)!=cadtemp.npos || cadtemp.find("+",0)!=cadtemp.npos)
				if (cadtemp[k+inc]=='-' || cadtemp[k+inc]=='+')
				{
					//Miro que hay a la izquierda del caracter, si hay un
					//espacio no pasa nada, sino, tengo que agregarle la e
					if (k>0) //Si no es el primer caracter
					{
						if (cadtemp[k-1]!=' ')
						{
							modifcoord[dim]=true;
							char temp=cadtemp[k];
							inc=inc+1;
							cadtemp[k]='E';
							cadtemp+=temp;
						}
					}
				}
			} //Delf for k

			//Leo las coordenadas del nodo 1 del CBUSH (o
			//que esta en la posicion del n`odo 2 o 3)
			vectemp[dim]=stringadouble(cadtemp);
		}//for dim

		//Ojo con el constructor de Nodos
		//construyo con id, id sistema de coordenadas, 
		Nodo n(id,idsc,vectemp);
		Nodos[indn]=n;
		//Ubicacion en NASTRAN
		//Nodos[indn].Nodo_Cargar_Ubic(i);
		indn++;
		//Nodos.push_back(n); //ESTO TARDA MUCHO!!!!;
	}// Recorrido de nodos

	return Nodos;
}

//Poner para varios tipos
//template <typename T>
//void Archivo::Leer_Campo(int campo, T&)
//{/
//
//}

//Leo sistemas de coordenadas
//Aca es una funcion de archivos que retorna un vector de sistema de coordenadas
const vector <SistCoord> Archivo::Leer_SistCoord()
{
	int i,j;

	vector <int> pos_sc;

	vector <SistCoord> Scs;


	////Chequeo donde empieza y terminan las filas de los SC
	pos_sc=Buscar_ini_fin("CORD");

	for (i=0;i<numfilas;i++)
	{
		if(lineas[i].find("CORD2",0)!=lineas[i].npos)
			filas_sistcoord.push_back(i);

	}
	int numsc=filas_sistcoord.size();

	string prueba;
	char cad[100];
	int tipo=1;
	int id;
	//Ahora tomo los datos de los sistemas de coordenadas
	for (i=0;i<numsc;i++)
	{
		sscanf(lineas[filas_sistcoord[i]].c_str(),"%s %d",cad,&id);
		prueba=cad;
		if (prueba=="CORD2R")
			tipo=0;
		if (prueba=="CORD2C")
			tipo=1;
		if (prueba=="CORD2S")
			tipo=2;

		SistCoord Sc(id,tipo);
		Scs.push_back(Sc);
	}

	return (const vector<SistCoord>) Scs;
}


int Archivo::Buscar_numnodos()
{
	int n=pos_nodos[1]-pos_nodos[0]+1;
	return n;
}

vector <string> Archivo::Cadena()
{
	return lineas;
}



void Archivo_sal::Escribir_Parte_cadena(const vector<string> cad,const vector <int> pos)
{
	for (int f=pos[0];f<pos[1];f++)
	{
		fsal<<cad[f]<<endl;
	}
}

void Archivo_sal::Escribir_cadenas(const vector<string> cad)
{
	int tam=cad.size();
	for (int f=0;f<tam;f++)
	{
		fsal<<cad[f]<<endl;
	}
}

void Archivo_sal::Escribir_cadena(const string cad)
{
	fsal<<cad<<endl;
}

Archivo_sal::Archivo_sal()
{}

void Archivo_sal::Iniciar(string cad)
{
	fsal.open (cad.c_str());
	//Abro archivo
	cout<< "Abriendo archivo de salida..."<<endl ;
	//fsal=fsalida;
}

const vector <int> Archivo::Pos_Nodos()	//Devuelve el inicio final de la posición de nodos
{
	return pos_nodos;
}

const int Archivo::MaxId_Mat()
{
	int maxid=0;
	int temp;
	for (int i=0;i<numfilas;i++)
	{
		if(lineas[i].find("MAT1",0)!=lineas[i].npos)
		{
			sscanf(lineas[i].c_str(),"%*s %d",&temp);
			if (temp>maxid)
				maxid=temp;
		}
	}
	return temp;
}

const int Archivo::MaxId_Entidad(string cad)
{
	int maxid=0;
	int temp;
	for (int i=0;i<numfilas;i++)
	{
		if(lineas[i].find(cad.c_str(),0)!=lineas[i].npos)
		{
			sscanf(lineas[i].c_str(),"%*s %d",&temp);
			if (temp>maxid)
				maxid=temp;
		}
	}
	return temp;
}


const int Archivo::NumFilas()
{
	return Archivo::numfilas;
}

void Archivo::Reemplazar_Linea(const int &pos, string cad)
{
	lineas[pos]=cad;
}

string Archivo::Linea(const int &pos)
{
	return lineas[pos];
}
