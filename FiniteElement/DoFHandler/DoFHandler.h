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

#ifndef _DOFHANDLER_H_
#define _DOFHANDLER_H_


#include "../Mesh/FEGrid.h"

#include "CSRGrid.h"	//TO MODIFY, WILL BE CONVERTED BY METIS
#include "metis.h"

namespace FluxSol
{

//DOF HANDLER
//AS IN DEAL.II, DOF HANDLER DIMENSION IS THE ELEMENT DIMENSION
//field_dim account for element degree of freemdom for each node, i.e. Element Degree
template <int dim>
class DoFHandler
{
    protected:
		int field_dim;		//Dim of the grid field (incognit)
        int dof_number;//AT FIRST IS NOT CONSTANT BECAUSE OF REMESH
		int node_number;

        const FeGrid <dim> & grid;
		std::vector < std::vector <unsigned int> > adjacent_dofs;
		std::vector <int> adj_dof_number;
		std::vector <unsigned int> nod_perm,nod_iperm;

		unsigned int max_row_width;

//        std::vector
//< types::global_dof_index > 	vertex_dofs
        //DIFFERENCE WITH NODAL DOF??
        std::vector <unsigned int> vertex_dofs;	//Like in deal.ii.
												//Here vertex_dofs values are simply the GLOBAL DOF
												//vector indices means NODE POSITIONS (not the id)
												//DOF Numbers is equal to mesh USEFUL node
		std::vector <unsigned int> node_pos;	//Inverse of vertex dofs, index are global dof, vector values
												//are node positions

		void Set_AdjDoFs();	//Internal


    public:
        //Constructors
        DoFHandler():grid(FeGrid<dim>()){}
        DoFHandler(const FeGrid<dim> &g);

        //TO MODIFY
        //SmartPointer< const
        //FiniteElement< dim, spacedim >
        //, DoFHandler< dim, spacedim > > 	selected_fe

        //TO MODIFY, WHICH OBJECT CALLS?
        //WHAT ABOUT INHERITED FROM ELEMENT?
        virtual void DistributeDoFs();

        inline const unsigned int & Vertex_DoF(const int &i)const{return this->vertex_dofs[i];}

    //Like in deal.ii
    inline const int  n_dofs(){return dof_number;};

        const std::vector <unsigned int> get_dof_indices (const Element<dim> &)const;
        const std::vector <unsigned int> get_dof_indices (const int &)const;

		const int & NumDoF() const {return this->dof_number;}

        //TO MODIFY
        //CAN BE A MEMBER FUNCTION OF A DIFFERENT CLASS
        //const vector <unsigned int> get_dof_indices (ITERATOR!!);
		const std::string outstr() const;

		const std::vector <int> & Adj_DoF_Number() const {return adj_dof_number;}

		const unsigned int & AdjDoF(const int &i, const int &j)const {return this->adjacent_dofs[i][j];}

		//inline const int & DOFPerm(const int &i){return (field_dim*nod_perm[)};


		const std::vector <unsigned int> RowWidth() const;

		void DoFReordering();

		const std::string RowWidth_outstr() const;
		void SetRowWidth();

		~DoFHandler(){}


};
//Constructor from grid
//At first global DOFs are asociated with Grid vertices
template <int dim>
DoFHandler<dim>::DoFHandler(const FeGrid<dim> &g)
:grid(g)
{
	cout << "Initializing DoF Handler..."<<endl;
	cout << "Num Elements: " <<g.NumElem()<<endl;
	cout << "Element field degree: " << endl;
	cout<<g.Elem(0).Degree();
	cout << "Field Dim:"<<endl;
	this->field_dim =g.Elem(0).Field_Dim();	//TO MODIFY, MUST BE OBTAINED BY THE ELEMENT
	cout << "Mesh degree: "  << g.Elref().Degree()<<endl;
    for (unsigned int n=0;n<g.NumNodes();n++) //TO MODIFY, SELECT ONLY USEFUL NODES!
	{
        this->vertex_dofs.push_back(n);
		//this->node_pos.push_back(n);
	}
	this->node_number=g.NumNodes();

	this->dof_number=vertex_dofs.size()*this->field_dim;
	cout << "vertex dof size" << vertex_dofs.size() << endl;
	cout << "field dim " << this->field_dim;
	cout << "DoF Number: " <<this->dof_number<<endl;
	cout << "adding adj dofs number" <<endl;
	std::vector<unsigned int> vt;
	adjacent_dofs.assign(this->dof_number,vt);
	//MUST TO LOCATE NONZERO ADJACENT DOFS FOR EACH ONE

	cout << "Element number " << g.NumElem() <<endl;
	for (unsigned int e=0;e<g.NumElem();e++)
	{
		//cout << "Creating element..." <<endl;
		Element <dim> el(g.Elem(e));
		//cout << "Element nodes" << el.NumNodes() <<endl;
		//cout << "Created." <<endl;
		for (int ne=0;ne<el.NumNodes();ne++)
		{
			for (int d1=0;d1<field_dim;d1++)	//Dimension
			{
				for (int ne2=0;ne2<el.NumNodes();ne2++)
				{
					for (int d2=0;d2<field_dim;d2++)	//Dimension
					{
						int dof=field_dim*el.NodePos(ne)+d1;
						int dof2=field_dim*el.NodePos(ne2)+d2;
						//cout << "dof 1 dof 2" << dof << " " <<dof2<<endl;
						adjacent_dofs[dof].push_back(dof2);
					}//d2
				}//ne2
			}//d1
		}//ne1
	}//element

	cout << "adding adj dofs number" <<endl;
	for (int dof=0;dof<this->dof_number;dof++)
	{
		adj_dof_number.push_back(adjacent_dofs[dof].size());

	}


	//Numbering
	//Sets renumbering
	CSRGrid <dim> csr (g);

	idx_t options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);

	// params_t *params=new params_t[1];

	  // params->ctype         = METIS_CTYPE_SHEM;
	  // params->iptype        = METIS_IPTYPE_NODE;
	  // params->rtype         = METIS_RTYPE_SEP1SIDED;

	  // params->ufactor       = OMETIS_DEFAULT_UFACTOR;
	  // params->pfactor       = 0;
	  // params->compress      = 1;
	  // params->ccorder       = 0;
	  // params->no2hop        = 0;

	  // params->nooutput      = 0;
	  // params->wgtflag       = 1;

	  // params->nseps         = 1;
	  // params->niter         = 10;

	  // params->seed          = -1;
	  // params->dbglvl        = 0;

	  // params->filename      = NULL;
	  // params->nparts        = 1;


	//According to ndmetis.
	options[METIS_OPTION_CTYPE]    = METIS_CTYPE_SHEM;

	options[METIS_OPTION_IPTYPE]   = METIS_IPTYPE_NODE;
	options[METIS_OPTION_RTYPE]    = METIS_IPTYPE_NODE;
	  options[METIS_OPTION_DBGLVL]   = 0;
	  options[METIS_OPTION_UFACTOR]  = 1;
	  options[METIS_OPTION_NO2HOP]   = 0;
	  options[METIS_OPTION_COMPRESS] = 1;
	  options[METIS_OPTION_CCORDER]  = 0;
	  options[METIS_OPTION_SEED]     = -1;
	  options[METIS_OPTION_NITER]    = 10;
	  options[METIS_OPTION_NSEPS]    = 1;
	  options[METIS_OPTION_PFACTOR]  = 0;

	  //TO MODIFY
	  options[16]  = 200;

	//for (int op=0;op<METIS_NOPTIONS;op++)
	//	cout << "METIS Option: " << op << " = " <<options[op]<<endl;

	idx_t *xadj=csr.XAdj();
	idx_t *adjncy=csr.Adjncy();

	//for (int n=0;n<g.NumNodes()+1;n++)
	//	cout << "csr xadj " << n << " = " <<xadj[n] <<endl;

	//for (int n=0;n<csr.AdjncySize();n++)
	//	cout << "adjncy " << " n " << adjncy[n]<<endl;


	idx_t *perm=new idx_t[g.NumNodes()];
	idx_t *iperm=new idx_t[g.NumNodes()];


	int error=METIS_NodeND((idx_t*)&g.NumNodes(), csr.XAdj(), csr.Adjncy(),NULL,options, perm, iperm);


	for (int n=0;n<g.NumNodes();n++)
	{
		nod_perm.push_back(perm[n]);
		nod_iperm.push_back(iperm[n]);
	}

    //for (int n=0;n<nod_perm.size();n++) cout << "perm iperm: " << nod_perm[n] << " " <<nod_iperm[n]<<endl;

	//BUT NO RENUMBERING

	void SetRowWidth();
	//for (int x=0;x<	g.NumNodes();x++) cout << "perm " << x << ": " << perm [x] <<endl;

	//for (int x=0;x<	g.NumNodes();x++) cout << "iperm" << x << ": " << iperm[x] <<endl;

	cout << "error : " <<error<<endl;


}


template <int dim>
void DoFHandler<dim>::DistributeDoFs()
{


    //deal.ii
  //if (dynamic_cast<const parallel::distributed::Triangulation<dim,spacedim>*>(&*tria) == 0)
  //  block_info_object.initialize(*this, false, true)

	for (int n=0;n<this->node_number;n++)
	{
		this->vertex_dofs[n]=this->nod_iperm[n];
	}

	this->Set_AdjDoFs();

}

//THIS CONVERTS
template <int dim>
const std::vector <unsigned int> DoFHandler<dim>::get_dof_indices (const int &e) const
{
    vector <unsigned int> v;

    for (int ne=0;ne<grid.Elem(e).NumNodes();ne++)
    {
        //tO MODIFY, WHERE MUST BE POINTING THIS INDICES, SINCE THE POSITION IS RELATIVE
        //
        for (int dir=0;dir<this->field_dim;dir++)
            v.push_back(field_dim*this->vertex_dofs[this->grid.Elem(e).NodePos(ne)] + dir);
    }

    const std::vector<unsigned int> cv(v);
    return cv;
    //return v;
}

template <int dim>
const std::string DoFHandler<dim>::outstr() const
{
	std::string str;
	std::ostringstream strs;

	str+="Dof Handler Global DOFs\n";
	strs.str("");
	for (unsigned int dof=0;dof<this->node_number;dof++)
	{
		strs<<this->vertex_dofs[dof]<<endl;
	}

	str+=strs.str();

	str+="Dof Handler Adjacent DOFs\n";
	for (unsigned int dof=0;dof<this->dof_number;dof++)
	{
		strs.str("");
		str+="Number of DOFS: ";
		strs<<" "<<std::setprecision(3)<<adjacent_dofs[dof].size()<<" ---- ";
		strs.str("");
		for (int dof2=0;dof2<adjacent_dofs[dof].size();dof2++)
			strs<<" "<<std::setprecision(3)<<adjacent_dofs[dof][dof2];
		str+=strs.str();
		str+="\n";
	}
	return str;
}

template <int dim>
void DoFHandler<dim>::Set_AdjDoFs()
{
	std::vector<unsigned int> vt;
	adjacent_dofs.clear();

	adjacent_dofs.assign(this->dof_number,vt);
	//MUST TO LOCATE NONZERO ADJACENT DOFS FOR EACH ONE
	for (unsigned int e=0;e<grid.NumElem();e++)
	{
		Element <dim> el(grid.Elem(e));
		for (int ne=0;ne<el.NumNodes();ne++)
		{
			for (int d1=0;d1<field_dim;d1++)	//Dimension
			{
				for (int ne2=0;ne2<el.NumNodes();ne2++)
				{
					for (int d2=0;d2<field_dim;d2++)	//Dimension
					{
						int dof =field_dim*this->vertex_dofs[el.NodePos(ne )]+d1;
						int dof2=field_dim*this->vertex_dofs[el.NodePos(ne2)]+d2;
						adjacent_dofs[dof].push_back(dof2);
						//cout << "Adj Dof Push [i][j]" << dof << " " <<dof2<<endl;
					}//d2
				}//ne2
			}//d1
		}//ne1
	}//element

	adj_dof_number.clear();
	for (int dof=0;dof<this->dof_number;dof++)
	{
		adj_dof_number.push_back(adjacent_dofs[dof].size());

	}


}


//At first RowWidth is not member
template <int dim>
const std::vector <unsigned int> DoFHandler<dim>::RowWidth() const
{
	std::vector<unsigned int > v;

	int row;
	for (unsigned int n=0;n<this->dof_number;n++)
	{
		unsigned int max=0;
		unsigned int min=1e8;

		for (int nn=0;nn<this->adj_dof_number[n];nn++)
		{
			if (adjacent_dofs[n][nn]>max)	max=adjacent_dofs[n][nn];
			if (adjacent_dofs[n][nn]<min)	min=adjacent_dofs[n][nn];
		}
		row = max-min+1;
		v.push_back(row);
	}

	return v;

}

template <int dim>
void DoFHandler<dim>::SetRowWidth()
{
	std::vector<unsigned int > v;
	max_row_width=0;
	int row;
	for (unsigned int n=0;n<this->dof_number;n++)
	{
		unsigned int max=0;
		unsigned int min=1e8;
		cout << "n: " << n <<endl;
		for (int nn=0;nn<this->adj_dof_number[n];nn++)
		{
			if (adjacent_dofs[n][nn]>max)	max=adjacent_dofs[n][nn];
			if (adjacent_dofs[n][nn]<min)	min=adjacent_dofs[n][nn];
		}
		row = max-min+1;
		if (row>max_row_width) max_row_width=row;
		v.push_back(row);
	}

	return v;

}

template <int dim>
const std::string DoFHandler<dim>::RowWidth_outstr() const
{
	std::vector <unsigned int> v =this->RowWidth();
	std::string str;
	std::ostringstream strs;

	str+="Row Width\n";


	for (unsigned int n=0;n<this->dof_number;n++)
		strs << "Row " << n << " : " <<v[n] << endl;

	str+=strs.str();

	unsigned int sum;

	strs.str("");
	strs<<"Max Row Width: " << this->max_row_width<<endl;


	for (unsigned int n=0;n<this->dof_number;n++) sum+=v[n];

	strs<<"Row Width Sum: " << sum << endl;
	str+=strs.str();

	return str;

}



} //FluxSol
#endif // _DOFHANDLER_H_
