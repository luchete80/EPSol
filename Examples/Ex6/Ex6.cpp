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
// Example 6//
//////////////

//Reordering DoF to reduce sparsity fill width
// Using METIS 

#include <petscksp.h>
#include "EPSol.h"

using namespace FluxSol;

int main(int argc,char **args)
{
	
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
	FluxSol::FeGrid<2> grid(1.0,1.0,0.,
							10,10,0);
							
	cout << "Total Grid Nodes" <<grid.NumNodes()<<endl;
	FluxSol::DoFHandler<2> dofhandler(grid);
	//cout << dofhandler.outstr();
	
	//MESH INFO
	cout << "Mesh Info" <<endl;
	cout << "Number of Elements" << grid.NumElem() <<endl;
	
	for (int e=0;e<grid.NumElem();e++)
	{
		cout << "Element nodes: ";
		for (int ne=0;ne<grid.Elem(e).NumNodes();ne++)
			cout << grid.Elem(e).NodePos(ne) << " " ;
		cout<< endl;
	}
	
	cout << "Initial AdjacentDoFs ..." << endl;
	
	cout << dofhandler.outstr();
	
	
	cout << dofhandler.RowWidth_outstr();
	
	cout << "Distributing DoFs..." << endl;
	
	dofhandler.DistributeDoFs();

	cout << dofhandler.RowWidth_outstr();
	cout << dofhandler.outstr();
	
	CSRGrid<2> csr(grid);
	csr.SavetoFile("out.graph");
	
	return 0;	
}