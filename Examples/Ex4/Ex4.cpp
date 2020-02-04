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

#include <string>
#include <iostream>
#include "SingleInputFile.h"
#include "EPSol.h"
#include "laspack.h"


int main()
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
	//c[0][0] = c[1][1] = ck;
	//c[0][1] = c[1][0] = ck*nu;
	//c[2][2] = ck*(1-nu)/2.;

	//Plain Strain
	ck = E*(1.-nu)/((1.+nu)*(1.-2*nu));
	c[0][0] = c[1][1] = ck;
	c[0][1] = c[1][0] = ck*nu/(1.-nu);
	c[2][2] = ck*(1 - 2*nu) / (2.*(1.-nu));

	//C = E / (1 - nu*nu)*[1 nu 0; nu 1 0; 0 0 (1 - nu) / 2];

	int totrows;
	Boolean Symmetry = False;

	size_t r;

	SetRTCAccuracy(1e-5);

	size_t cmp = 3;

	// INIT MESH

	//Mesh
	FluxSol::FeGrid<2> grid(input.Grid());
	cout << "Total Grid Nodes" <<grid.NumNodes()<<endl;
	FluxSol::DoFHandler<2> dofhandler(grid);
	
	//MESH INFO
	cout << "Mesh Info" <<endl;
	for (int e=0;e<grid.NumElem();e++)
		for (int n=0;n<4;n++)
			cout<<grid.Elem(e).NodePos(n)<<endl;

    //Global System
    QMatrix Kg;
	Vector Ug, Rg;

	totrows = dofhandler.NumDoF();
	cout << "Total DOFs" << totrows<<endl;
	Symmetry = False;
	Q_Constr(&Kg, "K", totrows, Symmetry, Rowws, Normal, True);
	V_Constr(&Ug, "U", totrows, Normal, True);
	V_Constr(&Rg, "R", totrows, Normal, True);

	//Dummy assembly Matrix
	FluxSol::Matrix<double> Kgi(dofhandler.NumDoF(),dofhandler.NumDoF());

    for (int r=1;r<=dofhandler.NumDoF();r++)
        Q_SetLen(&Kg, r, dofhandler.NumDoF());
		

	//Assemble of Matrix
	cout << "Assemblying matrix" <<endl;
	for (int e=0;e<grid.NumElem();e++)
	{
		cout <<"Element: " << e <<endl;
		// TO BE MODIFIED: DONT INCLUDE GRID
		//FluxSol::FEValues<2> fev(grid.Elem(e), grid);	
        vector <unsigned int > vn = dofhandler.get_dof_indices(e);
		
		for (int i=0;i<vn.size();i++)
			cout<<"GLOBAL DOFS"<<vn[i]<<endl;


		//FluxSol::FEValues<2> fev(qel,g);
		FluxSol::FEValues<2> fev(grid.Elem(e),grid);
		
        J = fev.Jacobian();
		B = fev.shape_grad_matrix();
		
		FluxSol::GaussMatrices H = grid.Elem(e).H();
		
		cout <<grid.XYZ(grid.Elem(e)).outstr();
		
		FluxSol::GaussFullMatrices dHdrs = fev.shape_localgrad_matrix();
		
		logfile << B.outstr();
		logfile << H.outstr();
		logfile << J.outstr();
		logfile << dHdrs.outstr();
		logfile << "Element Nodes Coordinates n\n";
		logfile << grid.XYZ(grid.Elem(e)).outstr();

        Kel.Clear();
		cout << "Creating Elemental Stiffness Matrix, GaussPoints: "<< intsch.NumPoints() <<endl;
        for (int g = 0; g < intsch.NumPoints(); g++)
        {

            FluxSol::Matrix<double> Kt = B.Mat(g).Tr()*c*B.Mat(g);
            for (int r = 0; r < 8; r++)
            for (int c = 0; c < 8; c++)
                Kel[r][c] += Kt[r][c]*intsch[g].w()*J.Mat(g).det();
			
			cout <<Kt.outstr();
        }
			cout <<"Element Matrix, r,c:" <<Kel.Rows()<<Kel.Cols()<<endl;
			cout << Kel.outstr()<<endl;
			cout <<"File output"<<endl;
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
                Kgi[vn[row]][vn[col]]+=Kel[row][col];
                //Q_SetEntry(&Kg, vn[col]+1, vn[row],vn[row]+1,Kel[row][col]);
            }

            //Assembly of load vector

        }

	}
	
	cout<<Kgi.outstr();

	//////////////////////////////////
	//Applying Boundary Condition ////
	//////////////////////////////////
    for (int row = 0; row < dofhandler.NumDoF(); row++)
        for (int col = 0; col < dofhandler.NumDoF(); col++)
            Q_SetEntry(&Kg, row+1, col,col+1,Kgi[row][col]);

	
	
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
				dof=2*node+d;
				for (int i = 0; i < dofhandler.NumDoF(); i++)
				{
					Q_SetEntry(&Kg, dof+1, i,i+1,0.0);
					Q_SetEntry(&Kg, i+1, dof,dof+1,0.0);
				}
				Q_SetEntry(&Kg, dof+1, dof,dof+1,1.0);				
			}
		}		
	
	}

    V_SetAllCmp(&Rg, 0.0);
	
	// node=2;
	// dir=0;
	// dof=2*node+dir;
		
		
	// V_SetCmp(&Rg, dof+1, -1000.0);
	
	
	for (int bcdof=0;bcdof<input.BLoadNumber();bcdof++)
	{
		//cout<<input.BLoad(bcdof).Show();

		node=input.BLoad(bcdof).NodeId();
		
		const vector<double>load=input.BLoad(bcdof).Load();
		//TO MODIFY, BETTER WITH A SINGLE VECTOR CONTAINING THE ONLY FIZED DIR
		for (int d=0;d<2;d++)
		{
				dof=2*node+d;
				V_SetCmp(&Rg, dof+1,load[d]);
		}		
	
	}


    for (int i = 0; i<20; i++)
	{
		V_SetAllCmp(&Ug, 0.0);
		BiCGSTABIter(&Kg, &Ug, &Rg, 1000, SSORPrecond, 1.0);
		cout << "It " << i << endl;
	}

	logfile << "\n Results \n";
	for (int j = 1; j < totrows+1; j++)
	{
		logfile << Ug.Cmp[j] << "\n";
	}

	Q_Destr(&Kg);V_Destr(&Rg);

	
	logfile.close();
	return 0;

}
