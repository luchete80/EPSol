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
#include "Modelo.h"
#include <vector>

using namespace std;
using namespace FluxSol::Nastran;
using namespace FluxSol;

////Constructor
Modelo::Modelo (string cad)
{
	int pos=cad.size()-4;
	string nombrelog=cad;
	nombrelog.replace(pos, 4, "-log");
	nombrelog+=".txt";
//	//
//	Nastran(cad);	//Construyo el archivo Nastran
	Nastran.Iniciar(cad);
	flog.open(nombrelog.c_str());
	//Inicio nodos
	numnodos=Nastran.Buscar_numnodos();
	//Inicio sistema de Coordenadas
	SistCoords=Nastran.Leer_SistCoord();
	numsistcoord=SistCoords.size();
	Nodos=Nastran.Leer_Nodos();
	Leer_Elementos();
	Asociar_Ids_Nodos();	//Esto tambien asocia las conectividades
	Asociar_Ids_Scs();

	inicio_mats_ec=false;
	inicio_prop_ec=false;
	inicio_carac_ec=false;

}

//Sólo busco posiciones de Nodos y elementos
void Modelo::Lectura_Rapida (string cad)
{
	int pos=cad.size()-4;
	string nombrelog=cad;
	nombrelog.replace(pos, 4, "-log");
	nombrelog+=".txt";
//	//
//	Nastran(cad);	//Construyo el archivo Nastran
	Nastran.Iniciar(cad);	//Esto ya busca la posicion de inicio y fin de los nodos
	flog.open(nombrelog.c_str());
	numnodos=Nastran.Buscar_numnodos();
	//Inicio nodos
}

Modelo::~Modelo()
{

}

void Modelo::Leer_Elementos()
{
	int elem=0;
	int icbush=0;
	int irigid=0;
	vector <int> lineas_cbush;
	vector <int> lineas_rigid;
	vector<int> conect;
	conect.assign(2,-1);

	cout <<"Leyendo Elementos"<<endl;
	for (int i=Nastran.pos_nodos[1]+1;i<Nastran.numfilas;i++)
	{
		//Elementos de fluidos
		if (Nastran.lineas[i].find("CQUAD4  ")!=Nastran.lineas[i].npos ||
			Nastran.lineas[i].find("CBEAM   ")!=Nastran.lineas[i].npos)
		{
			Elemento el;
			el.Leer_String(Nastran.lineas[i]);
			Elementos.push_back(el);
			elem++;
		}

	}//fin del for
	numelem=elem;

	//Asigno elementos
	int ind;
	int idprop;
	CBush temp;
	CBushes.assign(icbush,temp);	//Inicializo



	cout << "Asociando Indices de Elementos"<<endl;
	Asoc_Nodos.assign(MaxId_Nodos()+1,-1);
	for (int n=0;n<numnodos;n++)
	{
		Asoc_Nodos[Nodos[n].VerId_Nastran()]=n;
	}



	//Ahora actualizo los indices de elementos

}

//Asocia los ID's de NASTRAN con los internos
//Rehago los ids de acuerdo a los internos
//HAGO LAS DOS COSAS CREO UN INDICE DE ASOCIACIONES
void Modelo::Asociar_Ids_Nodos()
{
	cout << "Asociando Indices" <<endl;
	vector <int> v,conect;
	v.assign(2,-1);

	//Esto es muy lento
	//for (int i=0;i<numelem;i++)
	//{
	//	conect=Elementos[i].Conect();
	//	for (int n=0;n<numnodos;n++)
	//	{
	//		if (conect[0]==Nodos[n].VerId_Nastran())
	//		v[0]=n;
	//		else if (conect[1]==Nodos[n].VerId_Nastran())
	//		v[1]=n;
	//	}

	//	Elementos[i].Asignar_Conect_Interna(v);
	//}

	int nummax=0;
	int idn;
	//Ahora hago el indice grande
	for (int n=0;n<numnodos;n++)
	{
		idn=Nodos[n].VerId_Nastran();
		if (idn>nummax)
			nummax=idn;
	}

	//Le asigno el tamaño del máximo numero de NASTRAN +1
	IdNodoInterno.assign(nummax+1,-1);
	
	for (int n=0;n<numnodos;n++)
	{
		idn=Nodos[n].VerId_Nastran();
		IdNodoInterno[idn]=n;
	}

	for (int i=0;i<numelem;i++)
	{
		conect=Elementos[i].Conect();
		v[0]=IdNodoInterno[conect[0]];
		v[1]=IdNodoInterno[conect[1]];
		
		Elementos[i].Asignar_Conect_Interna(v);
	}
}

void Modelo::Asociar_Ids_Scs()
{
	cout << "Asociando Indices Sistemas de Coordenadas" <<endl;
	vector <int> v,conect;
	v.assign(2,-1);
	for (int i=0;i<numsistcoord;i++)
	{
		for (int n=0;n<numnodos;n++)
		{
			if (Nodos[n].Sc()==SistCoords[i].VerId())
			{
				Nodos[n].Cargar_IdSc_Int(i);
			}
		}
	}
}
// Inicio la longitud de los elementos TRUSS que van desde el nodo 2 a los 2x y 2y


void Modelo::Escribir_hasta_nodos()
{
	vector <string> lins=Nastran.Cadena();
	vector <int> pos=Nastran.Pos_Nodos();
	pos[0]=0;		//Quiero que escriba desde el principio del archivo
	pos[1]++;	//Falta el ultimo nodo
	Salida.Escribir_Parte_cadena(lins,pos);
}

void Modelo::Escribir_desde_fin_nodos()
{
	vector <string> lins=Nastran.Cadena();
	vector <int> pos=Nastran.Pos_Nodos();
	pos[0]=pos[1]+1;
	pos[1]=lins.size();
	Salida.Escribir_Parte_cadena(lins,pos);
}

//Requiere cargar el modelo
const int Modelo::MaxId_Nodos()
{
	int max=0;
	for (int n=0;n<numnodos;n++)
	{
		int id=Nodos[n].VerId_Nastran();
		if (max<id)
		max=id;
	}
	return max;
}

const int Modelo::MaxId_Nodos_Sin_Leer()
{
	int max=0;
	string cad;

	int id;
	int maxid=0;
	vector <int> pos=Nastran.Pos_Nodos();

	for (int i=pos[0];i<pos[1]+1;i++)
	{
		cad=Nastran.Linea(i);
		sscanf(cad.c_str(),"%*s %d",&id);
		if (id>maxid)	maxid=id;

	}
	return maxid;
}

const int Modelo::MaxId_Elem_Sin_Leer()
{
	int max=0;
	string cad;
	string tit;
	char tit2[100];
	int id;
	int maxid=0;
	vector <int> pos=Nastran.Pos_Nodos();

	for (int i=pos[1];i<Nastran.NumFilas();i++)
	{
		cad=Nastran.Linea(i);
		sscanf(cad.c_str(),"%s %d",tit2,&id);
		tit=tit2;
		if (tit.find("RBE2",0)!=tit.npos||tit.find("C",0)!=tit.npos)
		{
			if (id>maxid)	maxid=id;
		}
	}
	return maxid;
}

//Requiere cargar el modelo

const int Modelo::MaxId_Elems()
{
	int max=0;
	int id;
	for (int n=0;n<numcbush;n++)
	{
		id=CBushes[n].VerId();
		if (max<id)
		max=id;
	}

	for (int n=0;n<numrigid;n++)
	{
		id=Rigids[n].VerId();
		if (max<id)
		max=id;
	}

	return max;
}


void Modelo::Escribir_strings_Nodos_EdgeClamps()
{

	Salida.Escribir_cadenas(string_eclamps);
}

vector <string> Modelo::String_Eclamps()
{
	return string_eclamps;
}


vector <string> Modelo::String_PBUSH_Modif()
{
	return string_PBUSH_Modif;
}

vector <string> Modelo::String_CRods()
{
	return string_CRods;
}

vector <string> Modelo::String_RBEs()
{
	return string_RBEs;
}

vector <string> Modelo::String_RBEs_Centrales()
{
	return Modelo::string_RBEs_Centrales;
}

const int Modelo::Num_Cbush()
{
	return numcbush;
}

vector <string> Modelo::String_PRods()
{
	return string_PRods;
}

vector <string> Modelo::String_Mats()
{
	return string_Mats;
}
//Devuelve el id interno de Nastran
const int Modelo::VerId_CBush_Int(const int &idnastran)
{
	for (int cb=0;cb<numcbush;cb++)
	{
		if (CBushes[cb].VerId()==idnastran)
		return cb;
	}
}

//Leo la mu y Po de CBUSH
void Modelo::Leer_FriccionPrecarga_CBush(string cad)
{
	string s;
	ifstream fen(cad.c_str());
	double mu,po;
	int cb,idnas;
	int numcaract=0;
	bool enc;
	int fil=0;
	int idcarenc=0;
	bool encmu,encpo;
	while(getline(fen, s)) //
	{
		sscanf(s.c_str(),"%d %lf %lf",&idnas,&mu,&po);
		enc=true;
		for (int c=0;c<numcaract;c++)
		{
			//Tengo que ver si el mu y el po
			//enc significa que encontre un cbush con una caracteristica nueva
			if (mu==Caracts[c].VerMu() && po==Caracts[c].VerPo())
			{
				enc=false;
				idcarenc=c;
				break;
			}

		}
		if (fil==0) enc=true;

		if (enc)	//Nueva caract
		{
			idcarenc=numcaract;
			Caract Car(mu,po);
			Caracts.push_back(Car);
			numcaract++;
		}
		cb=VerId_CBush_Int(idnas);
		CBushes[cb].Iniciar_Caract(idcarenc);
		fil++;
	}

	numcaracts=numcaract;
	inicio_carac_ec=true;
}

//Sabiendo el vector Caract creo las propiedades de CBush
void Modelo::Iniciar_PRod_Segun_Caracteristicas(int &maxipbush)
{
	//Me fijo cual es el indice amyor de las propiedades
	if (!inicio_mats_ec)
	{
		cout << "Para Iniciar las Propiedades de los Edge Clamps deben iniciarse los materiales"<<endl;
		cout << "No se iniciaron Propiedades" <<endl;
		return;
	}
	Material m;


}

void Modelo::Iniciar_CRodsArea(const double &area)
{
	 crods_area=area;
}



void Modelo::Escribir_RBEs_Centrales()
{
	Salida.Escribir_cadenas(String_RBEs_Centrales());
}

void Modelo::Escribir_RBEs_Centrales_Alin_CRODs()
{
	Salida.Escribir_cadenas(String_RBEs_Alin_CRODS());
}

vector <string> Modelo::String_RBEs_Alin_CRODS()
{
	return string_RBEs_alin_CRODS;
}

void Modelo::Escribir_Materiales_Bilineal()
{
	for (int c=0;c<numcaracts;c++)
	{
		string linea=Mats[c].Cadena_Bilineal();
		Salida.Escribir_cadena(linea);
		Modelo::string_Mats.push_back(linea);
	}
}

void Modelo::Igualar_Agregados(const Modelo &m)
{

}

string Modelo::Linea_Nastran(const int &pos)
{
	return Nastran.Linea(pos);
}


void Modelo::Leer_Elementos_String()
{


    for (int i=Nastran.Pos_Nodos()[1];i<Nastran.NumFilas();i++)
    {

        Elemento el;

        string cad=Nastran.Linea(i);
        if ( (cad.find("CBUSH",0)!=cad.npos) || (cad.find("CQUADR",0)!=cad.npos) || (cad.find("CQUAD4",0)!=cad.npos) || (cad.find("CBEAM",0)!=cad.npos) ||
         (cad.find("CBAR",0)!=cad.npos) || (cad.find("CTRIAR",0)!=cad.npos) || (cad.find("CTRIA3",0)!=cad.npos) || (cad.find("RBE2",0)!=cad.npos) || (cad.find("RBE2",0)!=cad.npos)
            || (cad.find("CONM2",0)!=cad.npos))
        //Primero veo que elemento es
        {
        el.Leer_String(Nastran.Linea(i));
        Elementos.push_back(el);
        }
        //Si es un CBUSH pongo en el contenedor un derivado
//        else if (cad.find("CBUSH",0)!=cad.npos)
 //       {

   //     }

    }



}

//Lee todos los nodos de los rigid spider
void Modelo::Leer_Elementos_String_Mejorado()
{

    bool prim_linea;
    for (int i=Nastran.Pos_Nodos()[1];i<Nastran.NumFilas();i++)
    {
        vector <string> vcad;
        Elemento el;

        string cad=Nastran.Linea(i);
        vcad.push_back(cad);
        prim_linea=true;

        if ( (cad.find("CBUSH",0)!=cad.npos) || (cad.find("CQUADR",0)!=cad.npos) || (cad.find("CQUAD4",0)!=cad.npos) || (cad.find("CBEAM",0)!=cad.npos) ||
         (cad.find("CBAR",0)!=cad.npos) || (cad.find("CTRIAR",0)!=cad.npos) || (cad.find("CTRIA3",0)!=cad.npos) || (cad.find("RBE2",0)!=cad.npos) || (cad.find("RBE2",0)!=cad.npos)
            || (cad.find("CONM2",0)!=cad.npos))
        //Primero veo que elemento es
        {
            //Sumo las cadenas hasta que no aparezca el +
            //U otro elemento (cuya linea puede tener un +)
            //La primera comparacion es para la linea que tiene + al final, las que siguen lo tienen al ppio y al final
            int linea=0;
            while((prim_linea && (cad.find("+",0)!=cad.npos)) || ( (cad[0]=='+') && (!prim_linea)) )
            {
                prim_linea=false;
                i++;
                cad=Nastran.Linea(i);

                //Si no comenzo un nuevo elemento
                if (cad[0]=='+')
                {
                    vcad.push_back(cad);
                }
                else
                {
                i--;
                break;
                }
            }


        el.Leer_StringVec(vcad);
        Elementos.push_back(el);
        //cout <<i<<endl;
        }
        //Si es un CBUSH pongo en el contenedor un derivado
//        else if (cad.find("CBUSH",0)!=cad.npos)
 //       {

   //     }

    }



}

void Modelo::Buscar_Nodos_Aislados()
{
    //Creo un vector booleano con los indices reales de NASTRAN

    vector<Elemento>::iterator elit;
    vector<Elemento>::iterator elit2;
    bool encontro[2];

    //Recorro las conectividades
    //Podria recorrer solo los CBUSH
    int el=0;
    cout << "Buscando nodos sueltos " << el <<endl;
    for (elit=Elementos.begin();elit!=Elementos.end();elit++)
    {
        cout <<"Buscando elemento " <<el<<endl;
        bool suelto=true;

        vector<int> conect=(*elit).Conect();
        //Miro si los 2 nodos del CBUSH están unidos a algo
        if ((*elit).Es_cbush())
        {
            encontro[0]=encontro[1]=false; //Encontro nodo 1                //Para cada nodo del CBUSH recorro los elementos de nuevo y miro cada nodo
            //Repaso los 2 nodos del CBUSH
            for (int n=0;n<conect.size();n++)
            {


                for (elit2=Elementos.begin();elit2!=Elementos.end();elit2++)
                {
                    //Si el elemento no es el propio CBUSH
                    if (elit!=elit2)
                    {
                        vector<int> conect2=(*elit2).Conect();

                        //Si ya hay un nodo que

                        for (int n2=0;n2<conect2.size();n2++)
                        {
                            if (conect2[n2]==conect[n])
                            encontro[n]=true;
                        }
                    }

                }

                // Si ya no encontro el nodo 1 en ningun elemento no tiene sentido
                if (!encontro[0])
                    break;
            }//Fin de los dos nodos del CBUSH

            if (!encontro[0] || !encontro[1])
            {
                (*elit).Es_libre(true);
            }
            else{(*elit).Es_libre(false);}

        }//Fin de si es un CBUSH


        el++;
    }//Termino elit


}

//Convierto a Vertices de Modelo Fluido los nodos nastran
vector <_Vertex> & Modelo::Convert_VertexVector()
{
	vector <_Vertex> vv;

	for (int v=0;v<this->numnodos;v++)
	{
		_Vertex vert(this->Nodos[v].VerId_Nastran(),this->Nodos[v].Sc(),this->Nodos[v].Coords());
		vv.push_back(vert);
	}

	return vv;
}
