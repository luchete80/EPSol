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
#ifndef _FEGRID_H_
#define _FEGRID_H_

#include "Grid.h"
#include "Element.h"
#include "Node.h"

#include "Matrix\Matrix.h"

using namespace std;

namespace FluxSol{

// WHICH IS BEST??
// A TEMPLATE TYPE
// template <int dim, typename T < dim > >
// class pru
// {
	// std::vector  v <T>;
	// public:

// };

template<int dim>
class FeGrid:
public _Grid
{
	protected:
	unsigned int degree;
	//const Element<dim> &elref;
	vector <Node> node;
	//No const value, mesh can be modificated (remesh, ALE (inherited))
	//Elements are the cells
	std::vector< Element<dim> *> element;	//SEE INHERITANCE WITH CELL
	void Create (const double &lex, const double &ley, const double &lez,
		const int &nex, const int &ney, const int &nez);

	int num_nodes;

	public:
	//Constructors, MUST TO EXPLICIT ONE DIFFERET TO DEFAULT

		FeGrid(const std::vector<Node> &vn, const std::vector < Element<dim> > &ve)
		//:elref(ve[0])
		{
                this->num_cells = ve.size();
                cout << "ve size "<<num_cells<<endl;
                this->num_nodes = vn.size();

                //Including inherited object in base container
                for (int e = 0; e < this->num_cells; e++)
                    this->element.push_back(&ve[e]);

                for (int n = 0; n < this->num_nodes; n++)
                    this->node.push_back(vn[n]);

                cout << "Mesh Created " <<endl;
                cout << "Number of Nodes: " << this->NumNodes()<<endl;
                cout << "Number of Elements: " << this->num_cells <<endl;
		}


		FeGrid(const std::vector<Node> &vn, const std::vector < Element<dim> *> &ve)
		{

			this->num_cells = ve.size();
			cout << "ve size "<<num_cells<<endl;
			this->num_nodes = vn.size();

			//Including inherited object in base container
			for (int e = 0; e < this->num_cells; e++)
				this->element.push_back(ve[e]);

			for (int n = 0; n < this->num_nodes; n++)
				this->node.push_back(vn[n]);

		}

		FeGrid(const Element <dim> &el,
						const double &lex, const double &ley, const double &lez,
						const int &nex, const int &ney, const int &nez);

		FeGrid(const double &lex, const double &ley, const double &lez,
                const int &nex, const int &ney, const int &nez);


		FeGrid(){}
		~FeGrid(){}

		const Element<dim> & Elref()const{return *this->element[0];}

		FeGrid(const std::string &filename);

		inline const int & NumNodes(){return this->num_nodes;};

		//
		inline const Matrix<double> XYZ(const Element<dim> &el)const;


		inline const Node & Nod(const int &n)const{ return this->node[n]; }

        inline const Element <dim> & Elem (const int &e) const{return *this->element[e];}
        inline const int & NumElem(){return this->Num_Cells();}

		//inline operator=(const FeGrid<dim> right){}


};


	template<int dim>
	inline const Matrix<double> FeGrid<dim>::XYZ(const Element<dim> &el)const
	{
		Matrix<double> m(el.NumNodes(),dim);
		for (int d = 0; d < dim;d++)
		for (int n = 0; n < el.NumNodes(); n++)
		{
			Vec3D v = this->node[el.NodePos(n)].Coords();
			for (int i=0;i<3;i++)
			//cout << "vector coords"<< v[i] << " ";
			//cout <<endl;

			//Vec3D v = this->node[n].Coord[d];
			m[n][d] = v[d];
			//m[n][d] = this->node[n].Coords[d];

		}
		return m;
	}

	//NASTRAN FILE FORMAT
	template <int dim>
	inline FeGrid<dim>::FeGrid(const std::string &filename)
	//:elref(Element<dim>())
	{




	}


}//FluxSol


#endif
