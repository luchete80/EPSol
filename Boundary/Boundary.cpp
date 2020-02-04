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
#include "Boundary.h"

using namespace FluxSol;

//template <typename T,unsigned S> 
//Patch::Patch(T (&v)[S])
//{
//	//O sizeof [int]
//	int numfaces=arraysize(v);
//	for (int i=0;i<numfaces;i++)
//		this->id_face.push_back(facelist[i]);
//
//}


Patch::Patch(std::list <int> &idfaces)
{
	//O sizeof [int]
	std::list<int>::iterator it;
	this->num_faces=idfaces.size();
	for (it=idfaces.begin();it!=idfaces.end();it++)
		this->id_face.push_back((*it));

}

Boundary::Boundary (std::vector < Patch > patchlist)
{ 
	std::vector< Patch >::iterator it;
	num_faces=0;
	num_patches=0;
	for (it=patchlist.begin();it!=patchlist.end();it++)
	{
		patch.push_back(*it);
		this->num_faces+=it->Num_Faces();
		num_patches++;
	}

}