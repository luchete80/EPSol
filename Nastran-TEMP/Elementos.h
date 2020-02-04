#ifndef ELEMENTOS_H
#define ELEMENTOS_H

#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>


using namespace std;



class Elemento{
	public:
		Elemento();
		vector <int> Conect();
		vector <int> Conect_int();
		void Asignar_Conect_Interna(vector <int> v);
		void Asignar_Conect_Nastran(vector <int> v);
		const string Tipo(){return tipo;}

		void Leer_String(const std::string cad);
        void Leer_StringVec(const std::vector <string> cad);

		const int VerId();
		vector <int> VerConect_Int();
		Elemento* This();
		void Asignar_Prop(const int &pid);
		const int VerPid();
		void Elemento::Iniciar_Nastran(const int &i, vector <int> nodos_nastran);
		const int Elemento::Linea_Nastran();
		void Elemento::Iniciar_Fila_Nastran(const int &f);

		void Elemento::Iniciar_Prop(const int &i, const int &pos);
		void Elemento::Iniciar_Nodo_Master(const int &n);
		const int Elemento::Nodo_Master();
		const int Elemento::Pos_prop();
		int Elemento::Ver_Pos_Nastran();
		void Elemento::Iniciar_Prop_Modif(const int &v);
		void Elemento::Pos_Nastran(const int &p);
		const int Elemento::Ver_Prop_Modif();
        bool Es_cbush() {return es_cbush;}
        bool Es_libre() {return es_libre;}
        void Es_libre(bool es) {es_libre=es;}




		protected:	//Para que puedan acceder los elementos
		bool es_libre;
		bool es_cbush;
		int nodo_master;		//Indica cual nodo es el master
		int id;
		int nnodos;
		vector <int> Nodo;		//Indice del nodo NASTRAN
		vector <int> Nodos_int;	//Indice del nodo interno
		int pid;	//id de la propiedad, pero el interno
		int linea_nastran;	//Posicion de la fila en NASTRAN

		int id_nastran_prop_modificada;

		int idprop,pos_idprop_nastran;
		int pos_nastran;
        std::string tipo;

};


class CBush:public Elemento{

	bool es_estructural;
	public:
		CBush::CBush();
		CBush::CBush(const int &i, vector <int> nodos);
		void CBush::Iniciar_Es_Estructural(bool es);
		const bool CBush::Iniciar_Es_Estructural();
		void CBush::IniciarOri(int i);
		const int  CBush::VerOri();
		void CBush::Iniciar_Caract(const int &i);
		const int CBush::Ver_Caract();


	private:
	bool CBushCentral;	//Si está en el medio o en el extremo en la configuración
	int Ori;	//indice de la orientacion
	int Caract;	//Friccion y precarga
	int pos_nastran;
};

class Rigid:public Elemento{
	public:
		Rigid::Rigid();
		Rigid::Rigid(const int &i, vector <int> nodos,bool es_spider);
		bool Rigid::Es_Spider();
		string gdl;

	private:
	bool Spider;	//es simple o spider
};

#endif
