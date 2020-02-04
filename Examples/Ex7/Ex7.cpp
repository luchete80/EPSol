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
// SINGLE SQUARE 1 x 1 - PLAIN STRESS
// DISPLACEMENTS OF 0.001
// NONLINEAR EXAMPLE
// SMALL (AND LARGE??) DISPLACEMENTS
// I.E. TOTAL (AND UPDATED??) LAGRANGIAN FORMULATION

#include <petscksp.h>
#include "EPSol.h"

using namespace FluxSol;

int main(int argc,char **args)
{
	
	ofstream logfile;
	logfile.open("Logfile.txt");


	const FluxSol::FEIntegrationScheme intsch;

	FluxSol::GaussFullMatrices dHdrs;
	FluxSol::ShapeFunctionGroup shfngr;

	int dim=2;	

	FluxSol::Matrix<double> Kel(8, 8);
	FluxSol::Matrix<double> c(3, 3);
	


	double young = 206.0e9;
	double nu = 0.3;

	//Plain Stress
	double ck;
	
	ck = young / (1 - nu*nu);
	c[0][0] = c[1][1] = ck;
	c[0][1] = c[1][0] = ck*nu;
	c[2][2] = ck*(1-nu)/2.;

	//Plain Strain
	//ck = young*(1.-nu)/((1.+nu)*(1.-2*nu));
	//c[0][0] = c[1][1] = ck;
	//c[0][1] = c[1][0] = ck*nu/(1.-nu);
	//c[2][2] = ck*(1 - 2*nu) / (2.*(1.-nu));
	
	cout << "Material Matrix" << endl << c.outstr()<<endl;

	//C = E / (1 - nu*nu)*[1 nu 0; nu 1 0; 0 0 (1 - nu) / 2];

	int totrows;
	int row,col;

	// INIT MESH

	//Mesh
	cout << "Creating mesh..." <<endl;
	
	// TO MODIFY!! CORRECT
	//FluxSol::FeGrid<2> grid(1.0,1.0,0.,
	//						1,  1,  0);

	FluxSol::FeGrid<2> grid(Element<2>(1),
							1.0,1.0,0.,
							1,  1,  0);
	//MESH INFO
	cout << "Mesh Info" <<endl;
	cout << "Number of Elements" << grid.NumElem() <<endl;							
	cout << "Total Grid Nodes" <<grid.NumNodes()<<endl;
	
	FluxSol::DoFHandler<2> dofhandler(grid);
	//cout << dofhandler.outstr();
	
	
	for (int e=0;e<grid.NumElem();e++)
	{
		cout << "Element nodes: ";
		for (int ne=0;ne<grid.Elem(e).NumNodes();ne++)
			cout << grid.Elem(e).NodePos(ne) << " " ;
		cout<< endl;
	}
	
	cout << "Initial AdjacentDoFs ..." << endl;
	
	//cout << dofhandler.outstr();
	
	
	//cout << dofhandler.RowWidth_outstr();
	//NOT RENUMBERING FOR THIS SIMPLE EXAMPLE
	

	//Solver
	PETSC_Solver<double,2> solver(dofhandler.NumDoF());
	
	Matrix<double> KL(8,8);		//Linear element Matrix
	Matrix<double> KNL(8,8);		//Linear element Matrix
	
	
	// LOOK TROUGH LOAD STEPS
	int dtnum=10;
	double load, incload;
	load=1.03e8;
	incload=load/dtnum;
	load =0.;
	
	
	Vector<double> utot(dim*4);		//element TOTAL nodal displacements, u at time t
	Vector<double> uinc(dim*4);		//element INCREMENTAL nodal displacements
	
	
	//To Modify
	//Next FEField <dim> Utot
	Vector <double> Utot(dofhandler.NumDoF());	//Global accum displacement field, U at time t
	Vector <double> Uinc(dofhandler.NumDoF());	//Global accum displacement field 
	Vector <double> F(dofhandler.NumDoF());		//Global internal forces 
	
	uinc.Clear();
	
	cout << "Number of Time Increments: " << dtnum << endl;
	
	for (int dt=0;dt<dtnum;dt++)
	{
		cout << "Increment number " << dt <<endl;
		
		Utot+=Uinc;
		
		F.Clear();
		//Do not perform nodal redistribution
		cout << "Element number: " << grid.NumElem() <<endl;
		
		//TO MODIFY, NEXT MUST BE CALLED Solver.ClearMat
		for (int e=0;e<grid.NumElem();e++)
			for (int r = 0; r < 8; r++)
				for (int c = 0; c < 8; c++)
				{
					vector <unsigned int> vn = dofhandler.get_dof_indices(e);
					solver.SetMatVal(vn[r],vn[c],0.);
				}
		
		for (int e=0;e<grid.NumElem();e++)
		{
		
			//Pass Increment to elem
			
			utot+=uinc;
		
			vector <unsigned int> vn = dofhandler.get_dof_indices(e);
			cout << "Element DoF Indices "<<endl; for (int i=0;i<8;i++) cout <<vn[i] << " ";
			cout <<endl;
			
			
			for (int i=0;i<8;i++)	//TO MODIFY, NEXT RETURNED BY FeField
				uinc[i]=Uinc[vn[i]];

			for (int i=0;i<8;i++)	//TO MODIFY, NEXT RETURNED BY FeField
				utot[i]=Utot[vn[i]];
				
			FluxSol::FEValues<2> fev(grid.Elem(e), grid);
			FluxSol::GaussFullMatrices J = fev.Jacobian();
			FluxSol::GaussFullMatrices H = fev.shape_value_matrix();
		
			FluxSol::GaussFullMatrices BL0 = fev.shape_grad_matrix();
			
			

			Matrix<double> BL1(3,dim*grid.Elem(e).NumNodes());
			
			
			Matrix<double> BNL(2*dim,dim*grid.Elem(e).NumNodes());
			
			Vector<double> Fel(8);	//Elemental internal forces vector
			
			Matrix<double> BL (3,2*grid.Elem(e).NumNodes());
			

			
			FluxSol::GaussFullMatrices dHdrs = fev.shape_localgrad_matrix();

			FluxSol::ShapeFunctionGroup shfngr = grid.Elem(e).CreateShapeFunctionGroup();
			
			KL.Clear();
			Fel.Clear();
			
			Matrix<double> I(dim,dim); //Identity for Green Lagrange
			I.SetIdentity();

			//Calculate Deformation Gradient Tensor X
			//Calculate actual positions
			Matrix<double> xyz=grid.XYZ(grid.Elem(e));	//TO MODIFY, CHANGE TO ELEMENT FUNCTION
			for (int i=0;i<grid.Elem(e).NumNodes();i++)
				for (int d=0;d<2;d++)
					xyz[i][d]+=utot[2*i+d];
			
			cout << "xyz" <<endl;
			cout << xyz.outstr()<<endl;
					
			for (int g = 0; g < intsch.NumPoints(); g++)
			{
				FluxSol::Matrix<double> Hg = H.Mat(g);
				FluxSol::Matrix<double> l(2,2);
				
				for (int i=0;i<grid.Elem(e).NumNodes();i++)
				{
					//COULD BE H
					l[0][0]+=fev.shape_grad_component(i,g,0)*utot[2*i];
					l[0][1]+=fev.shape_grad_component(i,g,1)*utot[2*i];
					l[1][0]+=fev.shape_grad_component(i,g,0)*utot[2*i+1];
					l[1][1]+=fev.shape_grad_component(i,g,1)*utot[2*i+1];
				}
				
				cout << "utot"  << endl << utot.outstr()<<endl;
				cout << "Global Utot"  << endl << Utot.outstr()<<endl;
				cout << "Global Uinc"  << endl << Uinc.outstr()<<endl;
				
				for (int i=0;i<grid.Elem(e).NumNodes();i++)
				{
					
					BL1[0][2*i  ]=fev.shape_grad_component(i,g,0)*l[0][0];
					BL1[0][2*i+1]=fev.shape_grad_component(i,g,0)*l[1][0];
					
					BL1[1][2*i  ]=fev.shape_grad_component(i,g,1)*l[0][1];
					BL1[1][2*i+1]=fev.shape_grad_component(i,g,1)*l[1][1];
					
					BL1[2][2*i  ]=fev.shape_grad_component(i,g,1)*l[0][0]+fev.shape_grad_component(i,g,0)*l[0][0];
					BL1[2][2*i+1]=fev.shape_grad_component(i,g,1)*l[1][0]+fev.shape_grad_component(i,g,0)*l[1][1];
				}
				
				//cout << "Adding B" <<endl;
				BL=BL1+BL0.Mat(g);	//Or BL0[g]
				
				cout << "BL0 "<<endl;
				cout << BL0.Mat(g).outstr()<<endl;
				
				cout << "BL1 "<<endl;
				cout << BL1.outstr()<<endl;
				
				cout << "BL Gauss point "<< g <<endl;
				cout << BL.outstr();
				
				
				FluxSol::Matrix<double> Ki = BL.Tr()*c*BL;
				

				for (int r = 0; r < 8; r++)
					for (int c = 0; c < 8; c++)
						KL[r][c] += Ki[r][c]*intsch[g].w()*J.Mat(g).det();
						
				 cout << "gauss " << g <<endl;
				
					
				 Matrix<double> X = fev.shape_grad_comps().Mat(g)*xyz;
				 Matrix<double> E = (X.Tr()*X-I)*0.5;
				 
				 cout << "Green Tensor" <<endl;
				 cout << E.outstr();
				 
				 
				 
				 Vector<double> Ev(3);
				 Ev[0]=E[0][0];Ev[1]=E[1][1];Ev[2]=2*E[0][1];	//
				 
				 cout << "X tensor "<<endl;
				 cout << X.outstr();
				 
				 cout << "I tensor "<<endl;
				 cout << I.outstr();
				 
				// //Rearrangement of Green Lagrange
				
				Matrix<double> S (4,4);
				
				Vector<double> Sv= c*Ev;
				
				S[0][0]=Sv[0];S[1][1]=Sv[1];S[0][1]=Sv[2];S[1][0]=Sv[2];
				S[2][2]=Sv[0];S[3][3]=Sv[1];S[3][2]=Sv[2];S[2][3]=Sv[2];
				
				cout << "S Tensor" << S.outstr()<<endl;
				
				cout << "Sv Tensor" << Sv.outstr()<<endl;
				
				Vector<double> Fi=BL.Tr()*Sv;	//Internal Forces Integral 
				
				for (int r = 0; r < 8; r++)
					Fel[r]+=Fi[r]*intsch[g].w()*J.Mat(g).det();
					
				cout << "Fel" <<endl;
				cout << Fel.outstr() << endl;
				
				cout << "Fi" <<endl;
				cout << Fi.outstr() << endl;
				
				//NonLinear Deformation Matrix
				for (int i=0;i<grid.Elem(e).NumNodes();i++)
				{
					BNL[0][dim*i  ]=fev.shape_grad_component(i,g,0);
					BNL[1][dim*i  ]=fev.shape_grad_component(i,g,1);
					BNL[2][dim*i+1]=fev.shape_grad_component(i,g,0);
					BNL[3][dim*i+1]=fev.shape_grad_component(i,g,1);
				}
				
				Ki = BNL.Tr()*S*BNL;	

				for (int r = 0; r < 8; r++)
					for (int c = 0; c < 8; c++)
						KNL[r][c] += Ki[r][c]*intsch[g].w()*J.Mat(g).det();			
				
			}//End of integration points

			cout << "Linear Stiffness Matrix" <<endl<<KL.outstr();
			
			cout << "NonLinear Stiffness Matrix" <<endl<<KNL.outstr();

			Kel=KL+KNL;
			
			for (int c = 0; c < 8; c++)
				F[vn[c]]+=Fel[c];

			cout << "External Forces F" << endl;
			
			cout << F.outstr();
				
			cout <<endl;

			for (int r = 0; r < 8; r++)
				for (int c = 0; c < 8; c++)
					solver.AddMatVal(vn[r],vn[c],Kel[r][c]);
			
		}//Elements	
		
		/////////////////////////////////////////////
		//Applying Loads and Boundary Conditions ////
		/////////////////////////////////////////////
		solver.Flush();	

		
		// LOAD CONDITIONS		
		//RHS = R-F
		load+=incload;
		cout << "Load Value: " << load << endl;
		  
		for (int c = 0; c < dofhandler.NumDoF(); c++)
			//solver.SetbValues(c,0.);
			solver.SetbValues(c,-F[c]);
		
		//ADDING APPLIED FORCES
		//0 and 4 will be overritten by bc, however is a good practice to sum all vector
		solver.SetbValues(2,load-F[2]);
		solver.SetbValues(6,load-F[6]);
		     
	 	// BOUNDARY CONDITIONS
		solver.ApplyBCOnDoF(0, dofhandler);
		solver.ApplyBCOnDoF(1, dofhandler);
		solver.ApplyBCOnDoF(4, dofhandler);
		
		solver.SetbValues(3,0.);
		solver.SetbValues(5,0.);
		solver.SetbValues(7,0.);
		
		
		//solver.ApplyDispOnDoF(2, 0.001, dofhandler);
		//solver.ApplyDispOnDoF(6, 0.001, dofhandler);
	
		// SOLVE SYSTEM
		// NEXT TO MODIFY: LIKE OPENFOAM: Matrix K * FEField == FEField
		
		solver.ViewInfo();
		solver.Solve();
		solver.ViewInfo();
		
		//PASS SOLVER RESULTS TO FIELD
		//WILL BE NO NECESSARY
		Uinc=solver.X();
		
	
	}//End of dt
	
	return 0;	
}