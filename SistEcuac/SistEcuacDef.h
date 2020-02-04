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
#ifndef _SISTECUACDEF_H
#define _SISTECUACDEF_H


//In templates definitions must be included namespace instead of..
namespace FluxSol
{

template <typename T>
int FluxSol::Eqn<T>::Width()
{
	int min=1e10;
	int max=0;
	for (int i=0;i<this->neighbour_id.size();i++)
	{
		if (neighbour_id[i]<min)
			min=neighbour_id[i];
		if (neighbour_id[i]>max)
			max=neighbour_id[i];
	}
	if (this->id<min)	min=id;
	if (this->id>max)	max=id;
	return (max-min+1);
}

//Min Id of neighbour coefficient
template <typename T>
int FluxSol::Eqn<T>::MinNeigbourId()
{
	int min=1e10;
	for (int i=0;i<this->neighbour_id.size();i++)
	{
		if (neighbour_id[i]<min)
			min=neighbour_id[i];
	}
	if (this->id<min)	min=id;
	return min;
}

template <typename T>
void FluxSol::EqnSystem<T>::Insert (const FluxSol::Eqn <T> &ec)
{
	eqn.push_back(ec);
}



//Both ap must refer to the same id
//roght coeffs pass negative and left term pass as positive to rhs
//If neighbours refer to differents id, neighbour ids will be added
template <typename T>
Eqn<T> & Eqn<T>::operator ==(const Eqn<T> &right) const
{
	Eqn<T> ret;

	vector <int> tempright;

	tempright = right.neighbour_id;

	vector<int> ids;	//checked neighbours ids

	//Che if both central coeff are the same
	if (id == right.id)
	{
		//Must chek if both neighbour ids are the same and in that case
		ret.ap = this->ap - right.ap;

		bool lfound;
		//left eqn could have a bigger stencil 
		for (int nleft = 0; nleft < this->an.size(); nleft++)
		{
			lfound = false;
			int id;
			int nright = -1;
			while (!lfound && nright <= right.an.size())
			{
				nright++;
				if (this->neighbour_id[nleft] == right.neighbour_id[nright])
				{
					lfound = true;
					ids.push_back(nleft);
				}
			}

			//
			if (!lfound)
				ret.an[nleft] = this->an[nleft];
			else
			{
				ret.an[nleft] = this->an[nleft] - right.an[nright];
				tempright.erase(tempright.begin() + nright);
			}
				
		}//for nleft


		//Comparing right vectors
		bool rfound;
		for (int nright = 0; nright < tempright.size(); nright++)
		{
			//if ids[]
			rfound = false;
			int nleft = -1;
			while (!rfound)
			{
				nleft++;
				if (this->neighbour_id[nleft] == right.neighbour_id[nright])
				{
					rfound = true;
					ids.push_back(nleft);
				}
			}

			//if (!rfound)
			//	ret.an[ind] = this->an[ind];
			//else
			//	ret.an[ind] = this->an[ind] - right.an[nright];
		}

		//ret.An(this->Ap() - right.Ap());

		ret.source = right.source - this->source;

		return ret;
	}

}

//template <typename T>
//void EqnSystem<T>::Insert (const Eqn<T> &ec, std::vector<int> & nbreq)
//{
//	Insert(ec);
//	nbr_eqn=nbreq;
//}

// OPERATORS
template <typename T>
EqnSystem <T> & EqnSystem <T>::operator==(const double &d)
{
	Scalar source(d);

	//This value is added to each source term
	for (int e=0;e<eqn.size();e++)
	{
		Eqn(e)+=d;	//Overloaded Eqn operator, add to source term
	}

	return *this;
}


template <typename T>
void EqnSystem<T>::Log(std::string str)
{
	ofstream file(str);
	file<<"Equation System Info"<<endl;
	int e=0;
	for (int e=0; e<eqn.size();e++)
	{
		file<<"Eqn "<<e<<endl;
		eqn[e].Log(file);
	}
}

//template <typename T>
//EqnSystem<T>::EqnSystem(const Fv_CC_Grid &FvG)
//{
//	Fv_CC_Grid FvGrid;
//
//	dimension=int(pTraits<T>::rank);
//
//	for (int c=0;c<FvGrid.Num_Cells();c++)
//    {
//		vector <T> an;
//		T init_an;
//
//		//Insert empty Values
//		T ap,source;
//		vector <int> nbr_id;
//		an.assign(FvGrid.Cell(c).Num_Neighbours(),init_an),
//		nbr_id.assign(FvGrid.Cell(c).Num_Neighbours(),0);
//
//		for (int n=0;n<FvGrid.Cell(c).Num_Neighbours();n++)
//			nbr_id[n]=FvGrid.Cell(c).Neigbour(n);
//
//		FluxSol::Eqn <T> eq(c,ap,an,source,nbr_id);
//		eqn.push_back(eq);
//
//	}//End for cells
//
//}

//Binary members
//Arguments Can not be defined as constants
//Passing right coeffs to left and left source to sight
template <typename T>
EqnSystem <T> EqnSystem <T>::operator==(EqnSystem <T> &right)
{
	EqnSystem <T> ret;

	if (EqnV().size() == right.EqnV().size())
	{
		//Check if both eqns have same size
		for (int e = 0; e<EqnV().size(); e++)
		{
			ret.Eqn(e) = (eqn[e] == right.Eqn(e));
		}
	}

	return ret;
}


//
//template <typename T>
//EqnSystem<T>::EqnSystem(_CC_Fv_Field<T> &fvfield):field(fvfield)
//{
//	Fv_CC_Grid &FvGrid=this->field.Grid();
//
//	dimension=int(pTraits<T>::rank);
//
//	for (int c=0;c<FvGrid.Num_Cells();c++)
//    {
//		vector <T> an;
//		T init_an;
//
//		//Insert empty Values
//		T ap,source;
//		vector <int> nbr_id;
//		an.assign(FvGrid.Cell(c).Num_Neighbours(),init_an),
//		nbr_id.assign(FvGrid.Cell(c).Num_Neighbours(),0);
//
//		for (int n=0;n<FvGrid.Cell(c).Num_Neighbours();n++)
//			nbr_id[n]=FvGrid.Cell(c).Neigbour(n);
//
//		FluxSol::Eqn <T> eq(c,ap,an,source,nbr_id);
//		eqn.push_back(eq);
//
//	}//End for cells
//
//}


}

#endif