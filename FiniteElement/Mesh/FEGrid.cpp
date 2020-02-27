
#include "FEGrid.h"

using namespace FluxSol;

//TO MODIFY, TO SPECIFY ELEMENT TYPE IN GRID OR NOT
template <int dim>
FeGrid<dim>::FeGrid(const double &lex, const double &ley, const double &lez,
				const int &nex, const int &ney, const int &nez)
//:elref(Element<dim>())
{
	this->Create(lex,ley,lez,nex,ney,nez);
}

template <int dim>
void FeGrid<dim>::Create(const double &lex, const double &ley, const double &lez,
				const int &nex, const int &ney, const int &nez)
{
	double dx=lex/nex;double dy=ley/ney;
	this->num_nodes=(nex+1)*(ney+1);
	this->num_cells=0;

	double x,y;
	x=y=.0;
	vector <int> vn;	//Connectivities
	int n=0;
	for (int ny=0;ny<=ney;ny++)
	{
		x=0.;
		for (int nx=0;nx<=nex;nx++)
		{
			//cout << "Node"<<x<< " " << y <<endl;
			this->node.push_back(Node(n,x,y,0.0));
			x+=dx;
		}
		y+=dy;
	}

    //TO MODIFY, Z
	//cout << "nex, ney" << nex << " "<< ney<<endl;
	for (int ey=0;ey<ney;ey++)
	{
		for (int ex=0;ex<nex;ex++)		
		{
			//ifdef DEBUG
			//cout << "pushing back elems" <<endl;
			vn.clear();
			vn.push_back((nex+1)*(ey+1)+ex+1);
			vn.push_back((nex+1)*(ey+1)+ex);
			vn.push_back((nex+1)*ey + ex);
			vn.push_back((nex+1)*ey + ex+1);

			//Create an element from node vector
			//element.push_back(Element<dim>(vn));
			Element <dim> temp(this->degree,vn);
			Element <dim>* el=new Element <dim>;	//THIS IS NECESSARY IN ORDER TO ALLOCATE DATA, IF VECTOR IS FROM POINTERS
			*el=temp;
			//element.push_back(el);
			//TO MODIFY
			element.push_back(new Element <dim>(1,vn));

			this->num_cells++;
		}
	}
}

//TO MODIFY, CHANGE TO A POINTER SUCH AS ELEMENT TYPE CAN BE MODIFIED
template <int dim>
FeGrid<dim>::FeGrid(const Element <dim> &el,
				const double &lex, const double &ley, const double &lez,
				const int &nex, const int &ney, const int &nez)
//:elref(el)
{
	this->degree=el.Degree();
	this->Create(lex,ley,lez,nex,ney,nez);
	cout << "Mesh Created " <<endl;
	cout << "Number of Nodes: " << this->NumNodes()<<endl;
	cout << "Number of Elements: " << this->num_cells <<endl;
}



#include "FEGrid.inst"
