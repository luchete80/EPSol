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
#ifndef _CSRGRID_H_
#define _CSRGRID_H_

#include "FEGrid.h"
#include "metis.h"
#include <string>
#include <iostream>     // std::cout
#include <sstream>      // std::ostringstream
#include <iomanip>      // std::setprecision

namespace FluxSol
{
//Compress Storage Format Grid
template <int dim>
class CSRGrid
{
	//TO MODIFY: INERIT FROM PARENT BASIC GRID
	std::vector <idx_t> xadj;
	std::vector <idx_t> adjncy;

	idx_t *xadj2,*adjncy2;

	int num_nodes;
	int num_cells;

	public:
	//Constructors
	CSRGrid(const FeGrid<dim> &);
	const std::string outstr()const;
	idx_t * XAdj()const{return xadj2;};
	idx_t * Adjncy()const{return adjncy2;};

	const int & AdjncySize()const{return adjncy.size();};

	void SavetoFile(const std::string &) const;


};


//Function definitions
template <int dim>
CSRGrid<dim>::CSRGrid(const FeGrid<dim> &g)
{
	std::vector < std::vector <idx_t> > adj;
	//TO MODIFY, CHANGE TO STL map<> OR set <>
	std::vector<int> vv;
	adj.assign(g.NumNodes(),vv);
	this->num_nodes=g.NumNodes();
	this->num_cells=g.NumElem();
	cout  << "Creating csr" << endl;

	//cout << "Nodes " 	<< num_nodes <<endl;
	//cout << "Elements " << num_cells <<endl;

	//cout << "Creating adj...";
	//TO MODIFY
	//THIS DEPEND ON EACH ELEMENT, MUST HAVE TRAITS, AND MUST TO HAVE MAP
	for (unsigned int e=0;e<g.NumElem();e++)
	{
		Element <dim> el(g.Elem(e));

		//TO MODIFY: Grid function must to do this
		//TO MODIFY: See connectivity for other cell arrangements
		std::vector <int> elnodesext;
		elnodesext.assign(el.NumNodes()+2,0);
		for (int n=0;n<el.NumNodes();n++)
			elnodesext[n+1]=el.NodePos(n);

		elnodesext[0]=el.NodePos(el.NumNodes()-1);
		elnodesext[el.NumNodes()+1]=el.NodePos(0);

		//for (int n=0;n<el.NumNodes()+2;n++)
		//	cout << "Element Nodes Ext..." << elnodesext[n] <<endl;

		for (int ne=1;ne<el.NumNodes()+1;ne++)
		{
			int npos=elnodesext[ne];

			//cout << "Adding nodes.." << npos << ";"<<endl;

			bool found1=false;
			bool found2=false;
			for (int i=0;i<adj[npos].size();i++)
			{
				if (adj[npos][i]==elnodesext[ne+1]) found1=true;
				if (adj[npos][i]==elnodesext[ne-1]) found2=true;
			}


			if (!found1) adj[npos].push_back(elnodesext[ne+1]); /*cout <<"Adding adj "<<elnodesext[ne+1] << " ";*/
			if (!found2) adj[npos].push_back(elnodesext[ne-1]); /*cout <<elnodesext[ne-1] <<endl;*/
		}//ne1
	}//element

	//cout << "Creating xadj...";
	idx_t temp=0;
	xadj.push_back(temp);
	int pos=0;
	for (int n=0;n<g.NumNodes();n++)
	{
		for (int n2=0;n2<adj[n].size();n2++)
			adjncy.push_back(adj[n][n2]);

		temp=(idx_t)adj[n].size();
		pos+=temp;
		xadj.push_back(pos);
	}

	//cout << "Watching adj: "<<endl;
	//for (int i=0;i<num_nodes;i++)
	//{
	//	for (int j=0;j<adj[i].size();j++)
	//		cout << " " << adj[i][j];
	//	cout<<endl;
	//}

	xadj2=new idx_t[xadj.size()];
	adjncy2=new idx_t[adjncy.size()];

	for (int x=0;x<xadj.size();x++)
		xadj2[x]=xadj[x];

	for (int x=0;x<adjncy.size();x++)
		adjncy2[x]=adjncy[x];

//METIS_API(int) METIS_MeshToDual(idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind,
//                  idx_t *ncommon, idx_t *numflag, idx_t **r_xadj, idx_t **r_adjncy);

//METIS_API(int) METIS_MeshToNodal(idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind,
//                  idx_t *numflag, idx_t **r_xadj, idx_t **r_adjncy);


//int METIS NodeND(idx t *nvtxs, idx t *xadj, idx t *adjncy, idx t *vwgt, idx t *options,
//idx t *perm, idx t *iperm)

}

template <int dim>
const std::string
CSRGrid<dim>::outstr()const
{
	std::ostringstream strs;
	std::string str;

	strs.str("");
	strs<<"xadj Size: " <<xadj.size()<<endl;

	for (int n=0;n<this->num_nodes;n++)
		strs<< xadj[n] <<endl;

	str+=strs.str();
	strs.str("");

	strs<<"adjncy Size: " <<adjncy.size()<<endl;

	for (int n=0;n<this->adjncy.size();n++)
		strs<< adjncy[n] <<endl;

	str+=strs.str();

	return str;

}

template <int dim>
void CSRGrid<dim>::SavetoFile(const std::string &name) const
{
	ofstream outfile;
	outfile.open(name.c_str());

	outfile << this->num_nodes<<endl;

	for (int n=0;n<this->num_nodes;n++)
	{
		for (int nn=0;nn<(xadj[n+1]-xadj[n]);nn++)
			outfile<< this->adjncy[xadj[n]+nn] +1 <<" ";
		outfile <<endl;
	}
	outfile.close();
}

} //FluxSol
#endif
