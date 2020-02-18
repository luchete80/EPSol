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
// Plain Stress / Strain Problem
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
#include "laspack.h"

#include <iostream>
#include <vector>

using namespace std;

int main()
{

	ofstream logfile;
	logfile.open("Logfile.txt");

	//Element can be construted from vertex too
	//std::_Vertex

	std::vector<FluxSol::Node> v;

	//TO MODIFY: "ADDNODEFUNCTION IN GRID"
	v.push_back(FluxSol::Node(0, 1.0, 1.0, 0.0));
	v.push_back(FluxSol::Node(1, 0.0, 1.0, 0.0));
	v.push_back(FluxSol::Node(2, 0.0, 0.0, 0.0));
	v.push_back(FluxSol::Node(3, 1.0, 0.0, 0.0));

	const FluxSol::FEIntegrationScheme intsch;


	FluxSol::QuadLinearElement e(v);
	e.Set_Nodes(4, 0, 1, 2, 3);
	logfile << "Element Gauss Order: " << e.GaussOrder() << "\n\n";

	FluxSol::GaussMatrices H = e.H();

	std::vector< FluxSol::Element<2> > ve;
	ve.push_back(e);
	FluxSol::FeGrid <2> g(v, ve);
	FluxSol::FEValues<2> fev(e, g);

	FluxSol::GaussFullMatrices J = fev.Jacobian();

	FluxSol::GaussFullMatrices B = fev.shape_grad_matrix();

	FluxSol::GaussFullMatrices dHdrs = fev.shape_localgrad_matrix();

	FluxSol::ShapeFunctionGroup shfngr = e.CreateShapeFunctionGroup();

	double val = shfngr.ShapeFn(0).Val(0.577, 0.577, 0.);

	logfile << shfngr.ShapeFn(2).outstr();

	for (int d = 0; d < 4; d++)
	{
		logfile << "\n Shape Fn " << d << "\n";
		vector<FluxSol::ShapeFunction> df = shfngr.ShapeFn(d).diff();
		logfile << "\n Local Diff r Coeff \n";
		logfile << df[0].outstr();
		logfile << "\n Local Diff s Coeff \n";
		logfile << df[1].outstr();

	}

	logfile << "Jacobian Matrices \n\n";
	logfile << J.outstr();

	logfile << "Lineal Strain Matrices \n\n";
	logfile << B.outstr();

	logfile << "Shape Fn Matrices n\n";
	logfile << H.outstr();

	logfile << "Local grad Matrices \n\n";
	logfile << dHdrs.outstr();


	logfile << "Element Nodes Coordinates n\n";
	//logfile << g.XYZ(e).outstr();

	logfile << "Local derivative Functions n\n";
	logfile << fev.shape_grad_matrix().outstr();



	FluxSol::Matrix<double> Kel(8, 8);

	FluxSol::Matrix<double> c(3, 3);


	double E = 206.0e9;
	double nu = 0.3;

	//Plain Stress
	double ck = E / (1 - nu*nu);
	c[0][0] = c[1][1] = ck;
	c[0][1] = c[1][0] = ck*nu;
	c[2][2] = ck*(1 - nu) / 2.;

	//Plain Strain
	ck = E*(1. - nu) / ((1. + nu)*(1. - 2 * nu));
	c[0][0] = c[1][1] = ck;
	c[0][1] = c[1][0] = ck*nu / (1. - nu);
	c[2][2] = ck*(1 - 2 * nu) / (2.*(1. - nu));

	//C = E / (1 - nu*nu)*[1 nu 0; nu 1 0; 0 0 (1 - nu) / 2];

	for (int g = 0; g < intsch.NumPoints(); g++)
	{
		FluxSol::Matrix<double> Kg = B.Mat(g).Tr()*c*B.Mat(g);
		for (int r = 0; r < 8; r++)
		for (int c = 0; c < 8; c++)
			Kel[r][c] += Kg[r][c] * intsch[g].w()*J.Mat(g).det();

	}

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
	//////////////////// LASPACK SOLVER ///////////////////////


	QMatrix K;
	Vector U, R;

	int totrows = 8;
	Boolean Symmetry = False;
	Q_Constr(&K, "K", totrows, Symmetry, Rowws, Normal, True);
	V_Constr(&U, "U", totrows, Normal, True);
	V_Constr(&R, "R", totrows, Normal, True);

	size_t r;
	//Here Row begins from 1, but not column 
	for (r = 1; r < 9; r++)
	{
		Q_SetLen(&K, r, 8);
		for (size_t c = 0; c < 8; c++)
			//void Q_SetEntry(QMatrix *Q, size_t RoC, size_t Entry, size_t Pos, Real Val);
			//PArameter 3 is the sparse entry (from 0), pos is the real pos (from 1)
			Q_SetEntry(&K, r, c, c + 1, Kel[r - 1][c]);

	}


	V_SetAllCmp(&R, 0.0);
	V_SetAllCmp(&U, 0.0);
	SetRTCAccuracy(1e-5);

	size_t cmp = 3;
	V_SetCmp(&R, cmp, 1000.0);


	//SOLVER SYSTEM USING LASPACK
	for (int i = 0; i<20; i++)
	{
		V_SetAllCmp(&U, 0.0);
		BiCGSTABIter(&K, &U, &R, 1000, SSORPrecond, 1.0);
		cout << "It " << i << endl;
		//for (int j = 0; j<totrows; j++)
		//{
		//	Ui[j] = U.Cmp[j + 1];
		//	Ri[j] = R.Cmp[j + 1];
		//}

	}
	double a = U.Cmp[1];

	vector<double> Ui;
	Ui.assign(8, 0.);

	logfile << "\n Results \n";
	for (int j = 0; j < totrows; j++)
	{
		Ui[j] = U.Cmp[j + 1];
		logfile << Ui[j] << "\n";
	}


	Q_Destr(&K);


	logfile.close();
	return 0;


}