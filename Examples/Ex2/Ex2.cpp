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
// EXAMPLE 2 //
///////////////
#include "EPSol.h"

// Plain Stress / Strain Problem
// 10 x 10 Elements,4 node each, loaded as previous example
using namespace FluxSol;

int main()
{

	ofstream logfile;
	logfile.open("Logfile.txt");

	//Element can be construted from vertex too
	//std::_Vertex

	// std::vector<FluxSol::Node> v;

	// //TO MODIFY: "ADDNODEFUNCTION IN GRID"
	// v.push_back(FluxSol::Node(0, 1.0, 1.0,0.0));
	// v.push_back(FluxSol::Node(1, 0.0, 1.0, 0.0));
	// v.push_back(FluxSol::Node(2, 0.0, 0.0, 0.0));
	// v.push_back(FluxSol::Node(3, 1.0, 0.0, 0.0));

	const FluxSol::FEIntegrationScheme intsch(1,2);

	//Mesh
	FluxSol::FeGrid<2> grid(2.,1.,1.,2,1,1);
	cout << "Total Grid Nodes" <<grid.NumNodes()<<endl;

	FluxSol::DoFHandler<2> dofhandler(grid);

	// Element<2> e(grid.Elem(0));
	// e.Set_Nodes(4,0,1,2,3);
	// logfile << "Element Gauss Order: " << e.GaussOrder()<<"\n\n";

	// cout << "Assemblying matrix" <<endl;
	
	// FluxSol::GaussFullMatrices J,B,dHdrs;
	// FluxSol::ShapeFunctionGroup shfngr = e.CreateShapeFunctionGroup();

	// FluxSol::Matrix<double> Kel(8, 8);
	// FluxSol::Matrix<double> c(3, 3);

	// PETSC_Solver<double,2> solver(8);
    // Matrix<double> Kgi(dofhandler.NumDoF(),dofhandler.NumDoF());	

	// double E = 200.0e9;
	// double nu = 0.33;

	// //Plain Stress
	// double ck = E / (1 - nu*nu);

	// ck = E*(1.-nu)/((1.+nu)*(1.-2*nu));
	// c[0][0] = c[1][1] = ck;
	// c[0][1] = c[1][0] = ck*nu/(1.-nu);
	// c[2][2] = ck*(1 - 2*nu) / (2.*(1.-nu));
	

	// for (int e=0;e<grid.NumElem();e++)
	// {
		// cout <<"Element: " << e <<endl;
		// // TO BE MODIFIED: DONT INCLUDE GRID
		// //FluxSol::FEValues<2> fev(grid.Elem(e), grid);	
        // vector <unsigned int> vn = dofhandler.get_dof_indices(e);

		// for (int i=0;i<vn.size();i++)
			// cout<<"GLOBAL DOFS"<<vn[i]<<endl;
		
		// cout << "Creating values ..."<<endl;
		// FluxSol::FEValues<2> fev(grid.Elem(e),grid);
		
        // J = fev.Jacobian();
		// B = fev.shape_grad_matrix();
		
        // Kel.Clear();
		// cout << "Creating Elemental Stiffness Matrix, GaussPoints: "<< intsch.NumPoints() <<endl;
        // for (int g = 0; g < intsch.NumPoints(); g++)
        // {
			// cout << "B matrix " << B.Mat(g).outstr()<<endl;
			
            // FluxSol::Matrix<double> Kt = B.Mat(g).Tr()*c*B.Mat(g);
            // for (int r = 0; r < 8; r++)
				// for (int c = 0; c < 8; c++)
					// Kel[r][c] += Kt[r][c]*intsch[g].w()*J.Mat(g).det();
			
			// cout <<Kt.outstr();
        // }
		
		// logfile << "Stiffness Matrix \n\n";
		// logfile << Kel.outstr();

		// cout <<"Dof Indices"<<endl;
		// for(int row = 0; row < 8; row++)
			// cout <<vn[row]<<endl;
				
        // for (int row = 0; row < 8; row++)
        // {
            // for (int col = 0; col < 8; col++)
            // {
                // //Rows Assembly
                // //Laspack ,r from 1, c sparse matrix from 0, columntotal from 1
                // //Q_SetEntry(&Kg, r from 1, c from 0,c from 1,Kel[r][c]);
				// //PetscErrorCode  MatGetValues(Mat mat,PetscInt m,const PetscInt idxm[],PetscInt n,const PetscInt idxn[],PetscScalar v[])
					// //v	- a logically two-dimensional array for storing the values
					// //m, idxm	- the number of rows and their global indices
					// //n, idxn	- the number of columns and their global indices
				// //double val = 
				// //solver.MatVal(vn[row],vn[col],1);
                // Kgi[vn[row]][vn[col]]+=Kel[row][col];
				
				// solver.AddMatVal(vn[row],vn[col],Kel[row][col]);
            // }

        // }

	// }


	// logfile << shfngr.ShapeFn(2).outstr();

	// for (int d = 0; d < 4; d++)
	// {
		// logfile << "\n Shape Fn "<<d <<"\n";
		// vector<FluxSol::ShapeFunction> df = shfngr.ShapeFn(d).diff();
		// logfile << "\n Local Diff r Coeff \n";
		// logfile << df[0].outstr();
		// logfile << "\n Local Diff s Coeff \n";
		// logfile << df[1].outstr();

	// }






	// for (int g = 0; g < intsch.NumPoints(); g++){
		// FluxSol::Matrix<double> Kg = B.Mat(g).Tr()*c*B.Mat(g);
		// for (int r = 0; r < 8; r++)
		// for (int c = 0; c < 8; c++)
			// Kel[r][c] += Kg[r][c]*intsch[g].w()*J.Mat(g).det();}

	// logfile << "Stiffness Matrix \n\n";
	// logfile << Kel.outstr();


	// //Assemble of Matrix
	// for (int e=0;e<grid.NumElem();e++)
	// {
		
		// FluxSol::FEValues<2> fev(grid.Elem(e),grid);
        // vector <unsigned int > vn = dofhandler.get_dof_indices(e);

        // for (int i=0;i<8;i++)
            // cout<<vn[i]<<endl;

        // FluxSol::FEValues<2> feval(grid.Elem(e));

        // FluxSol::GaussFullMatrices J = fev.Jacobian();

        // FluxSol::GaussFullMatrices B = fev.shape_grad_matrix();

        // Kel.Clear();
        // for (int g = 0; g < intsch.NumPoints(); g++)
        // {
			// cout << "Gauss Point" << g<<endl;
            // cout << "Peso" << intsch[g].w();
            // FluxSol::Matrix<double> Kt = B.Mat(g).Tr()*c*B.Mat(g);
            // for (int r = 0; r < 8; r++)
            // for (int c = 0; c < 8; c++)
                // Kel[r][c] += Kt[r][c]*intsch[g].w()*J.Mat(g).det();

        // }

            // logfile << "Stiffness Matrix \n\n";
            // logfile << Kel.outstr();

        // for (int row = 0; row < 8; row++)
        // {
            // for (int col = 0; col < 8; col++)
            // {
                // //Rows Assembly
                // //Laspack ,r from 1, c sparse matrix from 0, columntotal from 1
                // //Q_SetEntry(&Kg, r from 1, c from 0,c from 1,Kel[r][c]);
                // Kgi[vn[row]][vn[col]]+=Kel[row][col];
                // //Q_SetEntry(&Kg, vn[col]+1, vn[row],vn[row]+1,Kel[row][col]);
            // }

            // //Assembly of load vector

        // }

	// }




	
	// logfile.close();
	// return 0;


}