#ifndef MODELO_H
#define MODELO_H

#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>

#include "./Nastran/Nastran.h"
#include "./Nastran/Nodo.h"
#include "./Nastran/Elementos.h"
#include "./Nastran/Caract.h"
#include "./Nastran/Prop.h"
#include "./Nastran/Material.h"

//Conversion a Volumenes Finitos
#include "./Boundary/Boundary.h"
//#include "./FiniteVolume/Mesh/FvGrid.h"

using namespace std;
//using namespace FluxSol;

namespace FluxSol{


namespace Nastran{
class Modelo{
	double long_rod;	//Long de los truss
	double crods_area;
	ofstream flog;

	int numnodos;
	int numelem;
	int numcbush;
	int numrigid;
	int numedgeclamps;
	vector <unsigned int> IdNodoInterno;	//vector de long de id nodo nastran maximo


	vector <Caract> Caracts;
	int numcaracts;

	int numsistcoord;

	//Vectores de elementos
	vector <SistCoord> SistCoords;



	vector <Material> Mats;


	//Boleanas para Edge Clamps
	bool inicio_mats_ec;
	bool inicio_prop_ec;
	bool inicio_carac_ec;


	bool bidim_model;

	void Modelo::Asociar_Ids_Nodos();



	//FUNCIONES PARA EDGECLAMPS

	//Las cadenas para la salida de edge clamps
	//Nodos y elementos
	vector <string> string_eclamps;	//SON LOS NODOS
	vector<string> string_PRods;
	vector<string> string_CRods;
	vector<string> string_RBEs;
	vector<string> string_Mats;

	//FUNCIONES
	void Modelo::Leer_Nodos();
	void Modelo::Leer_Elementos();


	public:

		Archivo Nastran;	//Aca ya le digo el constructor
		
		Modelo(){};
		Modelo(string cad);	//Lee de un archivo NASTRAN
		~Modelo();
		void Modelo::Lectura_Rapida (string cad);

		//Porque no a todos los modelos me interesa buscarle edge clamps, solo al reducido
		const int MaxId_Nodos_Sin_Leer();
		const int MaxId_Elem_Sin_Leer();
		void Buscar_Nodos_Aislados();


		const int Modelo::Num_Cbush();

		void Modelo::Leer_Coord_Rect(string  cad);
		void Modelo::Escribir_String_fsal(vector <string> cad);	//posiciones de filas iniciales y finales

		void Modelo::Modificar_Rigidez_Rigids();

		Archivo_sal Salida;	//Salida del modelo con sus edgeclamp

		void Modelo::Escribir_hasta_nodos();
		void Modelo::Asociar_Ids_Scs();
		const int Modelo::MaxId_Nodos();
		const int Modelo::MaxId_Elems();
		void Modelo::Escribir_strings_Nodos_EdgeClamps();
		void Modelo::CrearString_NodosEclamps2();
		void Modelo::Escribir_PBUSH_Modif();

		void Modelo::Escribir_desde_fin_nodos();

		//Creacion de Elementos
		void Modelo::Crear_CRods_EdgeClamps(int &maxid);
		void Modelo::Crear_Rigids_EdgeClamps(int &maxid);

		void Modelo::Leer_FriccionPrecarga_CBush(string cad);
		const int Modelo::VerId_CBush_Int(const int &idnastran);
		void Modelo::Iniciar_PRod_Segun_Caracteristicas(int &);

		//MATERIALES
		void Modelo::Iniciar_Mat_EdgeClamps(const int &maxmat, const double &huelgo, const double &po);
		void Modelo::Escribir_Materiales_Bilineal();

		void Modelo::Crear_RBEs_Centrales();
		void Modelo::Escribir_RBEs_Centrales();


		vector <string> Modelo::String_RBEs_Centrales();
		//Escritura en archivo
		void Modelo::Escribir_CRODS();
		void Modelo::Escribir_PRODS(const double &area);
		void Modelo::Escribir_RBEs_Eclamps();
		void Modelo::Crear_String_RBEs_Centrales(int &elem);

		void Modelo::Modificar_Rigidez_CBushes();
		void Modelo::Crear_String_RBEs_Centrales_a_CRODs(int &elem);	//Crea los rigids paralelos a los CRODS

		void Modelo::Igualar_Agregados(const Modelo &m);
		void Modelo::Escribir_RBEs_Centrales_Alin_CRODs();

		void Modelo::Reemplazar_Lineas(vector <string> cadenas);

		//DEVUELVO CADENAS PARA PASAR AL MODELO
		vector <string> Modelo::String_Eclamps();
		vector <string> Modelo::String_CRods();
		vector <string> Modelo::String_RBEs();
		vector <string> Modelo::String_PRods();
		vector <string> Modelo::String_Mats();
		vector <string> Modelo::String_RBEs_Alin_CRODS();
		vector <string> String_PBUSH_Modif();


		void Modelo::Listar_Caract_CBush();
		void Modelo::Iniciar_CRodsArea(const double &area);

		const int Modelo::Num_PBush();

		string Modelo::Linea_Nastran(const int &pos);
		void Modelo::Verificar_Dist_Nodo2_2xy(string cad, string sal);

		vector <CBush> CBushes;
		vector <Rigid> Rigids;	//Rigids Originales del Modelo
		vector <Rigid> Rigids_Eclamps;	//Rigids agregados
		vector <Rigid> Rigids_Eclamps_Centrales;
		vector <string> string_RBEs_Centrales;
		vector <string> string_RBEs_alin_CRODS;
		vector <string> string_PBUSH_Modif;	//nuevas propiedades de CBUSH

		void Modelo::Listar_Nodo_2_Insertos();
		vector <Nodo> Nodos_Eclamps;	//Lo hago en vector aparte

	void Modelo::Crear_String_PBUSH_Modif(int &maxipbush);
	void Modelo::Modificar_String_CBUSH();

	//Lo pongo publico para poder poner iteradores
	vector <Nodo> Nodos;
	vector <Elemento> Elementos;

	//Leo el string de ids y conectividades
	void Leer_Elementos_String();
    vector <int> Asoc_Nodos;	//la entrada es el ID de nastran y la salida el nodo
    void Leer_Elementos_String_Mejorado();  //Lee los slave de todos los nodos spring

	//Asociacion con Volumenes Finitos
	vector <_Vertex> & Convert_VertexVector();
	//void Modelo::Extract_Cells_and_Boundary(vector <Cell_CC> &ve, Boundary &bound);

	//Esta convierte la malla y llama a las anteriores
 	//Fv_CC_Grid & Convert2FvGrid();






};

}//Nastran
}//FluxSol

#endif
