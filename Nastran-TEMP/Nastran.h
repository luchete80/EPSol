#ifndef ARCHIVO_H
#define ARCHIVO_H

#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>

#include "./Mesh/Vertex.h"
#include "./Nastran/Cadenas.h"
#include "./Nastran/SistCoord.h"
#include "./Nastran/Nodo.h"
//#include "Modelo.h"

using namespace std;

namespace FluxSol {

//ARCHIVO DE ENTRADA!
class Archivo
{
        ifstream fid;		//Si el archivo es de entrada
        string nombre;
        ofstream fsal;	//Si el archivo es de salida


        vector <int> filas_sistcoord;	//Filas en las que está el string CORD2XXX

        //Funciones
        const vector<int> Buscar_ini_fin(std::string cadena);
        int Buscar_num_com(string cadena);
        void Buscar_pos_nodos();
        //Sistemas de coordenadas
		//	vector <SistCoord> Archivo::Leer_SistCoord();

        //vector <SistCoord> Archivo::Cargar_SistCoord(conect,ind_cbush,numcbush,numnodos);

	public:
		Archivo();
		//Archivo::Archivo(std::string cad);	//Constructor que abre el archivo

		void Iniciar(std::string cad);	//Esto es porque no se como iniciarlo
		int Buscar_numnodos();
		const vector <SistCoord> Leer_SistCoord();
		const vector <Nodo> Leer_Nodos();
		//const vector <Nodo> Archivo::Leer_Elementos();
		const vector <int> Pos_Nodos();	//Devuelve el inicio final de la posición de nodos

		//Esto es para probar que es mas rapido
		vector <string> lineas;
		vector <int> pos_nodos;	//Filas de inicio y fin
		int numfilas;
		void Escribir_Strings(vector <int> pos);	//posiciones de filas iniciales y finales

		vector <string> Cadena();

		//Con esto veo el maximo Id de material sin cargar los materiales
		const int MaxId_Mat();
		//Busca el maximo indice de cualquier entidad que este como PID en el campo 2
		const int MaxId_Entidad(string cad);

		void Reemplazar_Linea(const int &pos, string cad);
		string Linea(const int &pos);
		const int NumFilas();

     //   template <typename T>
     //   void Archivo::Leer_Campo(int campo, T&);




};

//ARCHIVO DE ENTRADA!
class Archivo_sal
{
	public:
		Archivo_sal();
		void Iniciar(string nombre);
		void Escribir_Parte_cadena(const vector<string> cad,const vector <int> pos);
		void Escribir_cadenas(const vector<string> cad);
		void Escribir_cadena(const string cad);
	private:

	ofstream fsal;	//Si el archivo es de salida



};

}

#endif
