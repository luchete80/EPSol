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
#ifndef _STRUCTURED_H
#define _STRUCTURED_H

#include "Grid.h"
#include <vector>

namespace FluxSol{
//Sirve para cualquier tipo de Malla, de volumenes finitos o de elementos finitos
//Con las variables ya iniciadas
//Se debio iniciar como una malla estructurada con sus numeros de vertices
//La plantilla es el tipo de malla
template<typename T>
void GenerarMallaRectangular(const double &lex,const double &ley, T& t)
{
	//Informacion de vertices
	int vx,vy;
	vector <int> verts_xyz=t.Num_Verts_xyz();
	vx=verts_xyz[0];
	vy=verts_xyz[1];

	int nex=verts_xyz[0]-1;
	int ney=verts_xyz[1]-1;


	double posx,posy,posz;

    posx=posy=posz=0.0;	
	//Creo vertices
    for (int k=0;k<2;k++)   //Vertices en z
    {
		posy=0.0;
        for (int j=0;j<vx;j++)	//Fila
        {
            for (int i=0;i<vy;i++)	//Columna
            {
                _Vertex pos(posx,posy,posz);				//Se va moviendo
                t.PushVert(pos);
                posx+=lex;
            }
            posy+=ley;
            posx=0.0;
        }
        posz+=1.0;
    }
//    //CASO Elemento_sS RECTANGULARES
//	//Se numeran como antes, de izquierda a derecha primero
//		//y luego de abajo arriba, NO ES COMO EN PAG 271
//
		int b=0;	//Son todos los vertices que ya pase al subir y
		vector <int> conect(8,-1);	//Vector de conectividades
		int num_verts_plano=vx*vy;
		for (int j=0;j<vy-1;j++)
		{
			for (int i=0;i<vx-1;i++)
			{
				//Aca los vertices se numeran como en pagina 352
				conect[0]=b+vx+i+1;
				conect[1]=b+vx+i;
				conect[2]=b+i;
				conect[3]=b+i+1;

				//Aca cambio de plano
				for (int k=0;k<4;k++)
				conect[k+4]=conect[k]+num_verts_plano;

                t.PushCell(conect);
			}//fin for de vx
			b+=vx;
		}	//fin for de vy

        //Debo generar el indice a los nodos, los indices de nodos son los mismos que los indices de cells
		t.Bidim_Mesh(true);
        t.Inicie_Cells(true);

}

}
#endif
