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

//////////////
// Example 5//
//////////////

//Solvig an elastic problem provided by input file, with petsc library
//Using PETSC_Solver class

#include <petscksp.h>
#include "EPSol.h"
#include "SingleInputFile.h"

using namespace FluxSol;

int main(int argc,char **args)
{

	cout <<"Opening file..."<<endl;
    SingleInputFile input("input.in");
	
	input.ShowData();
	
	ofstream logfile;
	logfile.open("Logfile.txt");


	const FluxSol::FEIntegrationScheme intsch;
	FluxSol::GaussMatrices H;

	FluxSol::GaussFullMatrices J,B,dHdrs;
	FluxSol::ShapeFunctionGroup shfngr;


	FluxSol::Matrix<double> Kel(8, 8);
	FluxSol::Matrix<double> c(3, 3);


	double E = 206.0e9;
	double nu = 0.3;

	//Plain Stress
	double ck;
	
	ck = E / (1 - nu*nu);
	c[0][0] = c[1][1] = ck;
	c[0][1] = c[1][0] = ck*nu;
	c[2][2] = ck*(1-nu)/2.;

	//Plain Strain
	//ck = E*(1.-nu)/((1.+nu)*(1.-2*nu));
	//c[0][0] = c[1][1] = ck;
	//c[0][1] = c[1][0] = ck*nu/(1.-nu);
	//c[2][2] = ck*(1 - 2*nu) / (2.*(1.-nu));

	//C = E / (1 - nu*nu)*[1 nu 0; nu 1 0; 0 0 (1 - nu) / 2];

	int totrows;
	int row,col;

	// INIT MESH

	//Mesh
	//To Modify, grid must be created with a reference element if must to be used with a dfhandler
	cout << "Input Grid Num Elem "<< input.Grid().NumElem()<<endl; 
	FluxSol::FeGrid<2> grid(input.Grid());
	cout << "Total Grid Nodes" <<grid.NumNodes()<<endl;
	
	cout << "Creating DoF Handler... " << endl;
	FluxSol::DoFHandler<2> dofhandler(grid);
	cout << dofhandler.outstr();
	
	//MESH INFO
	cout << "Mesh Info" <<endl;
	for (int e=0;e<grid.NumElem();e++)
		for (int n=0;n<4;n++)
			cout<<grid.Elem(e).NodePos(n)<<endl;

	//Dummy assembly Matrix
	FluxSol::Matrix<double> Kgi(dofhandler.NumDoF(),dofhandler.NumDoF());
	
	PETSC_Solver<double,2> solver(dofhandler.NumDoF());
	//TO MODIFY
	//Can be an integer or a non constant vector
	//solver.PreAllocateRows(dofhandler.Adj_DoF_Number());
	
	//Assemble of Matrix
	cout << "Assemblying matrix" <<endl;
	for (int e=0;e<grid.NumElem();e++)
	{
		cout <<"Element: " << e <<endl;
		// TO BE MODIFIED: DONT INCLUDE GRID
		//FluxSol::FEValues<2> fev(grid.Elem(e), grid);	
        vector <unsigned int> vn = dofhandler.get_dof_indices(e);

		for (int i=0;i<vn.size();i++)
			cout<<"GLOBAL DOFS"<<vn[i]<<endl;
		
		cout << "Creating values ..."<<endl;
		FluxSol::FEValues<2> fev(grid.Elem(e),grid);
		
        J = fev.Jacobian();
		B = fev.shape_grad_matrix();
		
        Kel.Clear();
		cout << "Creating Elemental Stiffness Matrix, GaussPoints: "<< intsch.NumPoints() <<endl;
        for (int g = 0; g < intsch.NumPoints(); g++)
        {
			cout << "B matrix " << B.Mat(g).outstr()<<endl;
			
            FluxSol::Matrix<double> Kt = B.Mat(g).Tr()*c*B.Mat(g);
            for (int r = 0; r < 8; r++)
				for (int c = 0; c < 8; c++)
					Kel[r][c] += Kt[r][c]*intsch[g].w()*J.Mat(g).det();
			
			cout <<Kt.outstr();
        }
		
		logfile << "Stiffness Matrix \n\n";
		logfile << Kel.outstr();

		cout <<"Dof Indices"<<endl;
		for(int row = 0; row < 8; row++)
			cout <<vn[row]<<endl;
				
        for (int row = 0; row < 8; row++)
        {
            for (int col = 0; col < 8; col++)
            {
                //Rows Assembly
                //Laspack ,r from 1, c sparse matrix from 0, columntotal from 1
                //Q_SetEntry(&Kg, r from 1, c from 0,c from 1,Kel[r][c]);
				//PetscErrorCode  MatGetValues(Mat mat,PetscInt m,const PetscInt idxm[],PetscInt n,const PetscInt idxn[],PetscScalar v[])
					//v	- a logically two-dimensional array for storing the values
					//m, idxm	- the number of rows and their global indices
					//n, idxn	- the number of columns and their global indices
				//double val = 
				//solver.MatVal(vn[row],vn[col],1);
                Kgi[vn[row]][vn[col]]+=Kel[row][col];
				
				solver.AddMatVal(vn[row],vn[col],Kel[row][col]);
            }

        }

	}
	
	cout<<Kgi.outstr();

	//////////////////////////////////
	//Applying Boundary Condition ////
	//////////////////////////////////
	solver.Flush();
	int node,dir,dof;
	for (int bcdof=0;bcdof<input.BFixNumber();bcdof++)
	{
		cout<<input.BFix(bcdof).Show();

		node=input.BFix(bcdof).NodeId();
		
		const vector<bool>fixed=input.BFix(bcdof).Fix();
		//TO MODIFY, BETTER WITH A SINGLE VECTOR CONTAINING THE ONLY FIZED DIR
		for (int d=0;d<2;d++)
		{
			if (fixed[d])
			{
				//If is a original nonzero value, must be added
				dof=2*node+d;
				for (int i = 0; i < dofhandler.Adj_DoF_Number()[dof]; i++)
				{
					int adjdof = dofhandler.AdjDoF(dof,i);
					
					//Have to set up only the nonzero matrix values
					solver.SetMatVal(dof,adjdof,0.);
					solver.SetMatVal(adjdof,dof,0.);
				}
				solver.SetMatVal(dof,dof,1.);		
			}
		}		
	
	}

	
	for (int bcdof=0;bcdof<input.BLoadNumber();bcdof++)
	{
		//cout<<input.BLoad(bcdof).Show();

		node=input.BLoad(bcdof).NodeId();
		
		const vector<double>load=input.BLoad(bcdof).Load();
		//TO MODIFY, BETTER WITH A SINGLE VECTOR CONTAINING THE ONLY FIZED DIR
		for (int d=0;d<2;d++)
		{
				dof=2*node+d;
				solver.SetbValues(dof,load[d]);
		}		
	
	}
	
	
	solver.ViewInfo();
	solver.Solve();
	solver.ViewInfo();

	return 0;	
}