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

///////////////
// EXAMPLE 1 //
///////////////
// Viwcoplastic Problem
// 1 Elements,4 node each, loaded as follows
//FORCEX
// ---------------
// |             |             
// |             |             
// |      1      |           
// |             |             
// |             |             
// ---------------
// FIX_XY		FIX_Y

//Solved with Laspack Library

#include "EPSol.h"

//Sparse libraries
//#include "laspack.h"

#include <iostream>
#include <vector>

using namespace std;
using namespace FluxSol;

class Eulerian_ViscoPlastic{

public:
    Eulerian_ViscoPlastic();
    void setup();
	void assemble();
	void calc_res();

private:
//Matrices
//Gauss
// ELEMENTAL MATRICES AND VECTORS
    GaussFullMatrices B,J,dHrs;
    
    vector< Matrix<double> > N;//Nv,Nsig,NF,Ns;
    Matrix<double> Bv,Bsig,BF,Bs,BL;

	FluxSol::Matrix<double> G,H,Kel;
	
	//Vectors
	//Solution
	Matrix<double> U,Uv,Usig,UF,Us;
	//Deformation gradients
	Matrix<double> F,P,L;	//Unsymmetric vectors
	
	//Residual 
	Matrix<double> R;	//Nv,Nsig,NF,Ns;
	
	//Dimensions
	int dim[4];
	enum field;
	
	//Material
	Matrix<double>c;
	double Ey,nu;
	

    ofstream logfile;
};


void Eulerian_ViscoPlastic::setup()
{

	ofstream logfile;
	logfile.open("Logfile.txt");
	
	//Matrices
	//Constants should be defined here
	G=H=Matrix<double>(4,4);
	G[0][0]=G[1][3]=G[2][3]=G[2][1]=1.;
	H[0][0]=G[0][1]=1.;
	
	Ey = 200.0e9;
	nu = 0.33;
	
	//Material
    c=Matrix<double>(4,4);
	double ck=Ey/((1.+nu)*(1.-2.*nu));
	c[0][0]=c[1][1]=c[2][2]=ck*(1.-nu);
	c[0][1]=c[0][2]=c[1][0]=c[2][0]=c[1][2]=c[2][1]=ck*(1.-nu);
	c[2][2]=ck*(0.5-nu);
	
	const FluxSol::FEIntegrationScheme intsch(1,2);


	//FluxSol::QuadLinearElement e(v);
	FluxSol::Element<2> e(v);
	e.Set_Nodes(4, 0, 1, 2, 3);
	//logfile << "Element Gauss Order: " << e.GaussOrder() << "\n\n";
	
	//Shape functions 
	//Nv = e.H();
	Bv=Matrix<double>(4,8);

	std::vector< FluxSol::Element<2> > ve;
	ve.push_back(e);
	FluxSol::FeGrid <2> g(v, ve);
    FluxSol::FEValues<2> fev(e, g);
    
	J = fev.Jacobian();
    B = fev.shape_grad_matrix();
   
    dHdrs=fev.shape_localgrad_matrix(); //COMPLETE
    ShapeFunctionGroup shfngr = e.CreateShapeFunctionGroup();
	// Shape Functions
	//To MODIFY REFER TO dim
    N[0]=Matrix<double>(2, 8);//v
	N[1]=Matrix<double>(4,16);//sig
	N[2]=Matrix<double>(4,16);//F
	N[3]=Matrix<double>(1, 4);//s

	//Solution (Fig 2.1)
	dim[0]=2;
	dim[1]=[2]=4;
	dim[3]=1;
	//R   =Matrix<double>(11,1);
	for (int i=0;i<4;i++)
		R[i]=Matrix<double>(dim[i],1);
	//Deformation gradients
	//Vector<double> F,P,L;	//Unsymmetric vectors
	
	//Residual 
	//Vector<double> Rv,Rsig,RF,Rs,R;	//Unsymmetric vectors
		


	Kel=Matrix<double> (44, 44);
}


void Eulerian_ViscoPlastic::assemble()
{

	for (int i=0;i<4;i++){
	//C = E / (1 - nu*nu)*[1 nu 0; nu 1 0; 0 0 (1 - nu) / 2];
	cout << "Num Integration Points"<<intsch.NumPoints()<<endl;
	for (int g = 0; g < intsch.NumPoints(); g++)
	{
		Bs=fev.shape_grad_comps(g);
		 
		Bv[0][2*i  ]=B[2][2*i]=Bs[0][i];
		Bv[1][2*i+1]=B[2][2*i]=Bs[0][i];
			
		FluxSol::Matrix<double> Kg = B.Mat(g).Tr()*c*B.Mat(g);
		for (int r = 0; r < 8; r++)
			for (int c = 0; c < 8; c++){
				Kel[r][c] += Kg[r][c] * intsch[g].w()*J.Mat(g).det();
				cout << "weight" << intsch[g].w()<<endl;
				cout << "Kel: "<< Kel[r][c]<<endl;	
	}}}


	logfile << "Stiffness Matrix \n\n";
	logfile << Kel.outstr();

	//DOFs 5, 6 and 8 restricted
	for (int i = 0; i < 8; i++)
	{
		Kel[4][i] = Kel[i][4] = 0.0;
		Kel[5][i] = Kel[i][5] = 0.0;
		Kel[7][i] = Kel[i][7] = 0.0;
	}
	Kel[4][4] = Kel[5][5] = Kel[7][7] = 1.;
	// //////////////////// LASPACK SOLVER ///////////////////////

	PETSC_Solver<double,2> solver(8);

	
	//DOFs 5, 6 and 8 restricted
	for (int i = 0; i < 8; i++)
	{
		Kel[4][i] = Kel[i][4] = 0.0;
		Kel[5][i] = Kel[i][5] = 0.0;
		Kel[7][i] = Kel[i][7] = 0.0;
	}
	Kel[4][4] = Kel[5][5] = Kel[7][7] = 1.;
	
	for (int row = 0; row < 8; row++){
		for (int col = 0; col < 8; col++)
		{
			solver.AddMatVal(row,col,Kel[row][col]);
		}

	}
	cout<<Kel.outstr();
	solver.Flush();
	
	solver.SetbValues(2,1000.0);
	

	solver.ViewInfo();
	solver.Solve();
	solver.ViewInfo();
	
	logfile.close();
	return 0;


}

void Eulerian_ViscoPlastic::calc_res()
{
	//R[0]=B[0].tr()*P - N[0].tr()*tP;
	// Rsig j = Bsig mjk * vk * Usig j - C mj E.(e)j   
	//R[1]=(N[1].tr()+tau(e)*v) * ()
	
}
