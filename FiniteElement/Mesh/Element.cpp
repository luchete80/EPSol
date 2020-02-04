#include "Element.h"

namespace FluxSol
{

//TO MODIFY: Make dregree inherited from poly space and constant
/* template <int dim>
Element<dim>::Element(const int &p)
{
	this->degree=p;
} */

//Can not return a reference because is not a preexistent variable
inline ShapeFunctionGroup QuadLinearElement::CreateShapeFunctionGroup() const
{
	ShapeFunctionGroup shfngr;

	//r,s,t,rs,st,rt,k
	//CAN MAKE A GENERAL FUNCTION BUT IN ONLY THIS GROUP DOES NOT OCCUPY SPACE
	shfngr.AddShapeFn(ShapeFunction( 0.25, 0.25, 0., 0.25, 0., 0., 0.25));//(1+r)(1+s)/4
	shfngr.AddShapeFn(ShapeFunction(-0.25, 0.25, 0.,-0.25, 0., 0., 0.25));//(1-r)(1+s)/4
	shfngr.AddShapeFn(ShapeFunction(-0.25,-0.25, 0., 0.25, 0., 0., 0.25));//(1-r)(1-s)/4
	shfngr.AddShapeFn(ShapeFunction( 0.25,-0.25, 0.,-0.25, 0., 0., 0.25));//(1+r)(1-s)/4


	return shfngr;

}


//ShapeFunctionGroup Q4TH::CreateShapeFunctionGroup() const
//{
//	ShapeFunctionGroup shfngr;
//
//	//r,s,t,rs,st,rt,k
//
//    shfngr.AddShapeFn(ShapeFunction( 0.25, 0.25, 0., 0.25, 0., 0., 0.25));//(1+r)(1+s)/4
//    shfngr.AddShapeFn(ShapeFunction(-0.25, 0.25, 0.,-0.25, 0., 0., 0.25));//(1-r)(1+s)/4
//    shfngr.AddShapeFn(ShapeFunction(-0.25,-0.25, 0., 0.25, 0., 0., 0.25));//(1-r)(1-s)/4
//    shfngr.AddShapeFn(ShapeFunction( 0.25,-0.25, 0.,-0.25, 0., 0., 0.25));//(1+r)(1-s)/4
//
//	return shfngr;
//}

//inline ShapeFunctionGroup H8TH::CreateShapeFunctionGroup() const
//{
//	ShapeFunctionGroup shfngr;
//
//    //Coeff Order
//	//r,s,t,rs,st,rt,k, rst
//    shfngr.AddShapeFn(ShapeFunction( 0.125, 0.125,-0.125, 0.125,-0.125,-0.125, 0.125,-0.125));//(1+r)(1+s)(1-t)/8
//    shfngr.AddShapeFn(ShapeFunction(-0.125, 0.125,-0.125,-0.125,-0.125, 0.125, 0.125, 0.125));//(1-r)(1+s)/4
//    shfngr.AddShapeFn(ShapeFunction(-0.125,-0.125,-0.125, 0.125, 0.125, 0.125, 0.125,-0.125));//(1-r)(1-s)/4
//    shfngr.AddShapeFn(ShapeFunction( 0.125,-0.125,-0.125,-0.125, 0.125,-0.125, 0.125, 0.125));//(1+r)(1-s)/4
//    shfngr.AddShapeFn(ShapeFunction( 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125));//(1+r)(1-s)/4
//    shfngr.AddShapeFn(ShapeFunction(-0.125, 0.125, 0.125,-0.125, 0.125,-0.125, 0.125,-0.125));//(1+r)(1-s)/4
//    shfngr.AddShapeFn(ShapeFunction(-0.125,-0.125, 0.125, 0.125,-0.125,-0.125, 0.125, 0.125));//(1+r)(1-s)/4
//    shfngr.AddShapeFn(ShapeFunction( 0.125,-0.125, 0.125,-0.125,-0.125, 0.125, 0.125,-0.125));//(1+r)(1-s)/4
//
//	return shfngr;
//}

//TO BE MODIFIED
//THIS MUST TO USE TRAITS
template <int dim>
inline ShapeFunctionGroup Element<dim>::CreateShapeFunctionGroup() const
{
	ShapeFunctionGroup shfngr;

	//r,s,t,rs,st,rt,k
	//CAN MAKE A GENERAL FUNCTION BUT IN ONLY THIS GROUP DOES NOT OCCUPY SPACE
	if (dim==2)
	{
		shfngr.AddShapeFn(ShapeFunction( 0.25, 0.25, 0., 0.25, 0., 0., 0.25));//(1+r)(1+s)/4
		shfngr.AddShapeFn(ShapeFunction(-0.25, 0.25, 0.,-0.25, 0., 0., 0.25));//(1-r)(1+s)/4
		shfngr.AddShapeFn(ShapeFunction(-0.25,-0.25, 0., 0.25, 0., 0., 0.25));//(1-r)(1-s)/4
		shfngr.AddShapeFn(ShapeFunction( 0.25,-0.25, 0.,-0.25, 0., 0., 0.25));//(1+r)(1-s)/4
	}

	return shfngr;
}


//INLINE??
template <int dim>
GaussMatrices Element<dim>::H() const	//Calls
{
	//	//Initializes ShapeFunction Matrix
	ShapeFunctionGroup shfngr = this->CreateShapeFunctionGroup();

	GaussIntegrationScheme intsch(this->GaussOrder(),dim);
	//	//For each gauss point create each sparse matrix
	GaussMatrices ret(dim, this->DOF(), dim*shfngr.Size(), intsch);

	Matrix<double> temp;

	for (int g = 0; g < intsch.NumPoints(); g++)
	{
		for (int h = 0; h < shfngr.Size(); h++)
		{
			//This MUST BE A MEMBER FUNCTION
			double fh = shfngr[h].Val(intsch[g]);
			ret[g][h] = fh;
			ret[g][h + shfngr.Size()] = fh;
			//Update locations
		}
	}

	for (int f = 0; f < shfngr.Size(); f++)
	{
		ret.SetPair(f, 0, 2 * f); ret.SetPair(f + shfngr.Size(), 1, (2 * f) + 1);
	}
	//	//Now i have element function values
	//
	//	//depends on nodes degrees of freedom and type
	//	//QUE SEA VIRTUAL
	//	//VIRTUAL e.ActiveShapeDOF() returning a sparsecoeff

	//SPECIALIZATION
	//THIS CAN BE GENERIC

	return ret;
}

template <int dim>
GaussMatrices Element<dim>::dHdrst() const
{
	//NON VIRTUAL
	//GENERIC (FOR ALL INHERITED) ELEMENTS

	//GaussMatrices ret(3, shfngr.Size(), intsch);
	//GENERAL
	ShapeFunctionGroup shfngr = this->CreateShapeFunctionGroup();

	GaussIntegrationScheme intsch(this->GaussOrder(),dim);

	GaussMatrices ret(dim, this->DOF(), dim*shfngr.Size(), intsch);

	vector<ShapeFunction> vsh;

	//PARTICULAR - VIRTUAL - PART
	for (int h = 0; h < shfngr.Size(); h++)
	{
		//vsh = shfngr.ShapeFn(h).diff();
		for (int g = 0; g < intsch.NumPoints(); g++)
		{
			//vsh=shfngr.
			//This MUST BE A MEMBER FUNCTION
			double fh = shfngr[h].Val(intsch[g]);
			ret[g][h] = fh;
			ret[g][h + shfngr.Size()] = fh;
			//Update locations
		}
	}


	return ret;

}


//IF NOTSPECIFIED
template<int dim>
Element<dim>::Element(const std::vector<Node> &vn)
:gaussor(quad),
field_dim(dim)
	//:gaussorder(2)
{
	this->degree=dim;
	int gaussor = linear;

	this->num_nodes = vn.size();

	//Assign indexes
	for (int n = 0; n < this->num_nodes; n++)
	{
		id_node.push_back(vn[n].Id());

	}

}




//RETURN NODE COORDINATES
template <int dim>
const Matrix<double> Element<dim>::XYZ()const
{

	Matrix<double> ret(this->num_nodes,dim);

    //return ret;
}

template<int dim>
void Element<dim>::Set_Nodes(const int numnodes, ...)
{

	va_list vl;
	va_start(vl, numnodes);
	//v_.push_back(va_arg(vl, T));
	for (int i = 0; i < numnodes; ++i)
	{

			this->id_node.push_back(va_arg(vl, int));

	}
	va_end(vl);

}

template <int dim>
Element<dim>::Element(const std::vector<int> &v)
:field_dim(dim)//default
{
	this->degree=(unsigned int) dim;	//TO MODIFY
	for (int i=0;i<v.size();i++)
		this->id_node.push_back(v[i]);
	this->num_nodes=v.size();
}

template <int dim>
void Element<dim>::AddConnectivity(const std::vector<int> &v)
{
	for (int i=0;i<v.size();i++)
		this->id_node.push_back(v[i]);
	this->num_nodes=v.size();
}

template <int dim>
Element<dim>::Element(const int &p, const std::vector<int> &v)
:Element(p)
{
	this->AddConnectivity(v);
}


}

#include "Element.inst"
