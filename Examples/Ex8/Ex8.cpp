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
// Example 8//
//////////////

// THERMAL EXAMPLE
// 1 x 1 Dim Square
// Single example with Zero temperature at bottom left corner
// And Heat Flux of 1 at Upper right corner


#include <petscksp.h>
#include "EPSol.h"

using namespace FluxSol;

int main(int argc,char **args)
{
	
	ofstream logfile;
	logfile.open("Logfile.txt");


	const FEIntegrationScheme intsch;

	FluxSol::GaussFullMatrices dHdrs;
	FluxSol::ShapeFunctionGroup shfngr;

	int dim=2;	

	FluxSol::Matrix<double> Kel(4, 4);
	

	int totrows;
	int row,col; 

	// INIT MESH

	//Mesh
	cout << "Creating mesh..." <<endl;
	FluxSol::FeGrid<2> grid(Element<2>(1),
							1.0,1.0,0.,
							2,  2,  0);

	//MESH INFO
	cout << "Mesh Info" <<endl;
	cout << "Grid Dim " <<  grid.Elref().Degree()<<endl;
	cout << "Number of Elements" << grid.NumElem() <<endl;							
	cout << "Total Grid Nodes" <<grid.NumNodes()<<endl;
		
	for (int e=0;e<grid.NumElem();e++)
	{
		cout << "Element " << e << "nodes: "<<endl;
		for (int ne=0;ne<grid.Elem(e).NumNodes();ne++)
			cout << grid.Elem(e).NodePos(ne) << " " ;
		cout<< endl;
	}
	
	cout << "Creating DoF Handler..."<<endl;
	FluxSol::DoFHandler<2> dofhandler(grid);
	//cout << dofhandler.outstr();
	

	
	cout << "Initial AdjacentDoFs ..." << endl;
	
	cout << dofhandler.outstr();
	
	
	cout << dofhandler.RowWidth_outstr();
	//NOT RENUMBERING FOR THIS SIMPLE EXAMPLE
	

	//Solver
	PETSC_Solver<double,2> solver(dofhandler.NumDoF());


	//Do not perform nodal redistribution
	cout << "Element number: " << grid.NumElem() <<endl;

	
	for (int e=0;e<grid.NumElem();e++)
	{
		vector <unsigned int> vn = dofhandler.get_dof_indices(e);
		cout << "vn size: " << vn.size() <<endl;
		cout << "Element DoF Indices "<<endl; for (int i=0;i<4;i++) cout <<vn[i] << " ";
		cout <<endl;

		FluxSol::FEValues<2> fev(grid.Elem(e), grid);
		FluxSol::GaussFullMatrices J = fev.Jacobian();
		FluxSol::GaussFullMatrices H = fev.shape_value_matrix();
	
		FluxSol::GaussFullMatrices B = fev.shape_grad_matrix();

		FluxSol::GaussFullMatrices dHdrs = fev.shape_localgrad_matrix();

		FluxSol::ShapeFunctionGroup shfngr = grid.Elem(e).CreateShapeFunctionGroup();
		
		Kel.Clear();
		cout << "Creating Elemental Stiffness Matrix, GaussPoints: "<< intsch.NumPoints() <<endl;
        for (int g = 0; g < intsch.NumPoints(); g++)
        {
			cout << "B Matrix " << B.Mat(g).outstr()<<endl;
            FluxSol::Matrix<double> Kt = B.Mat(g).Tr()*B.Mat(g);
            for (int r = 0; r < 4; r++)
            for (int c = 0; c < 4; c++)
                Kel[r][c] += Kt[r][c]*intsch[g].w()*J.Mat(g).det();
        }
		
		
		cout << "Elemental Matrix" << Kel.outstr();

		for (int r = 0; r < 4; r++)
			for (int c = 0; c < 4; c++)
				solver.AddMatVal(vn[r],vn[c],Kel[r][c]);
		
	}//Elements	
	
	solver.ViewInfo();
	
	/////////////////////////////////////////////
	//Applying Loads and Boundary Conditions ////
	/////////////////////////////////////////////
	solver.Flush();	
	  
	//ADDING APPLIED FORCES
	//0 and 4 will be overritten by bc, however is a good practice to sum all vector
	solver.SetbValues(dofhandler.NumDoF()-1,1.);
		 
	// BOUNDARY CONDITIONS
	solver.ApplyBCOnDoF(0, dofhandler);
		
	
	//solver.ApplyDispOnDoF(2, 0.001, dofhandler);
	//solver.ApplyDispOnDoF(6, 0.001, dofhandler);

	// SOLVE SYSTEM
	// NEXT TO MODIFY: LIKE OPENFOAM: Matrix K * FEField == FEField
	
	solver.ViewInfo();
	solver.Solve();
	solver.ViewInfo();
		
	
	return 0;	
}