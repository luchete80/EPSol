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
#define TRACEON
#include "Grid.h"
//#include "Log.h"

using namespace std;
using namespace FluxSol;
//Recordar que no son lo mismo que nodos
void _Grid::Num_Verts(const int &nv)
{
	num_verts=nv;
}

//Funcion que crea un grid de clase BASE con los vertices.
//Todo lo que es nodos se define en la clase derivada
_Grid::_Grid(const int &nex, const int &ney)
{
        int nv;
        num_cells_xyz.push_back(nex  );num_cells_xyz.push_back(ney  );num_cells_xyz.push_back(1);
        num_verts_xyz.push_back(nex+1);num_verts_xyz.push_back(ney+1);num_verts_xyz.push_back(2);
        nv=(nex+1)*(nex+2)*2;
        Num_Verts(nv);
        Num_Verts_xyz(nex,ney,2);
		Num_Cells(num_cells_xyz[0]*num_cells_xyz[1]*num_cells_xyz[2]);
}

//Construct from cells
_Grid::_Grid(std::vector<_Cell> vc)
{
	this->num_cells = vc.size();

	for (int c = 0; c < this->num_cells; c++)
	{
		cell.push_back(vc[c]);
	}

}

//Pasa los argumentos geometricos al constructor del Face
//Devuelvo vertices
vector <_Vertex> & _Grid::Verts(vector<int> indvert)
{
    vector <_Vertex> verts; //Vertices que componen la cara
    _Vertex sm;             //Punto medio de la superficie
    _Vertex sn;             //Vector superficie (incluye la dimension de la cara)
    //Recorro los indices
    //Se setean la posicion media de
    for (int i=0;i<indvert.size();i++)
    {
        //Voy llenando con los vertices que me vienen dando
        verts.push_back(vert[i]);
    }
    return verts;
    //Llamo al constructor de face
}


//INICIO CARAS
//ESTA ES UNA FUNCION MUY LARGA, DEBERIA PARTIRLA
//Esta funcion inicia las caras a partir de una malla que tiene vertices y conectividades
void _Grid::Iniciar_Caras()
{
//	int v1,v2;
//
//	int nf=0;
//	int nfi=0;	//Numero de caras interiores
//	int nfb=0;	//Numero de caras en el borde
//	bool enc;
//	int dc=0;		//dummy cells
//	int c2;
//
//    //Declaro otro iterador que no esta aca
//
//    //Inicio las celdas
    //if (inicie_cells)
    //{
    //    //Esto antes era sin iteradores
    //    //Ahora recorro asi las celdas
    //    for (cellit=cell.begin();cellit!=cell.end();cellit++)
    //    {
    //        //Esto lo puedo poner en otra funcion
    //        for (int i=0;i<4;i++)//primer vertice local
    //        {
    //            //Tomo los 2 primeros vertices
    //            v1=(*cellit).Id_Vert(i);
    //            if (i<3)	v2=(*cellit).Id_Vert(i+1);
    //            else		v2=(*cellit).Id_Vert(0);

    //            enc=false;
    //            vector <_Cell>::iterator cellit2;
    //            for (cellit2=cell.begin();cellit2!=cell.end();cellit2++)//Recorrida por todos los elementos para ver quien comparte lado
    //            {
    //                for (int i21=0;i21<4;i21++)
    //                    if ((*cellit2).Id_Vert(i21)==v2)//Encontre 2 elementos con 1 vertice comun
    //                        //Cuando los indices de un elemento gira hacia un lado
    //                        //los del vecino giran en sentido opuesto
    //                        //Busco otro vertice comun
    //                    if ((*cellit2).Id_Vert(i21+1)==v1 || (i21==3 && (*cellit2).Id_Vert(0)==v1))//Encontre el otro lado
    //                                                                                     //Esta comparacion
    //                                                                                //descarta a un cell consigo mismo
    //                    {
    //                        enc=true;
    //                        //Siempre el numero de cara se identifica CON EL PRIMER VERTICE
    //                        if(c2>(*cellit).Id()) //Poruqe de lo contrario esta cara ya esta contabilizada
    //                        {
    //                            (*cellit2).face[i21]=nf;
    //                            (*cellit).face[i]=nf;	//Es la cara que corresponde al primer vertice

    //                            //Aca cambio esto y construyo la cara
    //                            Face[nf].BF=false;
    //                            Face[nf].Cell[0]=c;		//ESTO ESTA HECHO DE FORMA TAL QUE EL ELEMENTO 0 SEA
    //                            Face[nf].Cell[1]=c2;	//DESDE DONDE SALE EL VECTOR DE LAS S
    //                            Face[nf].Vert[0]=v1;Face[nf].Vert[1]=v2;

    //                            nfi++;nf++;
    //                        }
    //                    }
    //            }// for de c2

    //            if (!enc)	//La cara es del borde
    //            {
    //                //Cell[num_cells_int+dc].Dummy=true;
    //                face[nf].Cell[1]=num_cells_int+dc;	//SIEMPRE LA CELDA 1 ES LA DE DUMMY!!!!!!
    //                face[nf].Cell_int[0]=c; Face[nf].Cell_ext[0]=num_cells_int+dc;
    //                face[nf].Cell[0]=c;
    //                //Falta numerar la variable Face del Cell dummy!!!!
    //                //El cell ext es c2
    //                (*cellit).Face[i]=nf;	//Es la cara que corresponde al primer vertice
    //                Face[nf].BF=true;
    //                Face[nf].Vert[0]=v1;Face[nf].Vert[1]=v2;
    //                nfb++;nf++;
    //                dc++;
    //            }

    //        }//Fin del for del vertice del elemento 1

    //        //PENDIENTE ACA AUMENTAR A DOS CELDAS

    //    }//Fin for cells

    //}//Fin de si ya inicie celdas
//
//	num_faces_int=nfi;
//	num_faces_ext=nfb;
//	num_faces_int=nf;
//
//	//ATENCION! esto es si un dummy
//	num_cells_d=num_faces_ext*dummyxface;	//los dummy son los numero de caras
//	num_cells=num_cells_d+num_cells_int;
//
//	//Iniciadas las ceras inicio los indices de stencils para los vertices
//	//LO TENGO QUE HACER EN UN FOR APARTE PORQUE TENGO QUE INICIAR TODOS LOS CELLS Y DE LA OTRA FORMA
//	//podria haber hecho un vector temporal pero bueh
//	int cv;
//	for (int c=0;c<num_cells_int;c++)
//	{
//		cell[c].id=c;      //Antes era esto
//		for (int f=1;f<5;f++)//recorro caras
//		{
//			if (face[Cell[c].face[f-1]].Cell[0]==c)
//				cell[c].Cell[f]=face[Cell[c].Face[f-1]].Cell[1];
//			else
//				cell[c].Cell[f]=face[Cell[c].Face[f-1]].Cell[0];
//		}
//	}
}


// LOG DE ARCHIVO
void _Grid::Log_Verts() // Imprime log de Vertices
{
    //recorro los vertices con el iterador perteneciente a la clase
    for (vertit=vert.begin();vertit!=vert.end();vertit++)
    {
        //meshlog<<(*vertit).Imprimir_Coord()<<endl;
    }

}


//El virtual solo se coloca en la declaracion!
void _Grid::Log()
{
    //Ojo que el log anda mal
//    ofstream flog("Hola.txt");
//    Logvecf(num_cells_xyz,flog);
//    Tracevec(num_verts_xyz);
//    LOG(num_verts)

}

void _Grid::Trace()
{
//    Tracevec(num_cells_xyz);
//    Tracevec(num_verts_xyz);
//    LOG(num_verts)
}
