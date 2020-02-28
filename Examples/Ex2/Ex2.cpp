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

namespace FluxSol{
class outfile
{
	protected:
	//To check if new fields have same mesh
	//Fv_CC_Grid &grid;
	FluxSol::FeGrid<2> &grid;
	vector <Vec3D> node_data;   //Results
	//vector <> cell_data;		//Element results

	public:
	//Constructors
	outfile(string name, FluxSol::FeGrid<2> &,const vector<Vec3D> &);

};
}

int main()
{

	ofstream logfile;
	logfile.open("Logfile.txt");
	

	FEIntegrationScheme intsch(1,2);
	

	//Mesh
	//FeGrid(const double &lex, const double &ley, const double &lez,
	int nex=100;
	int ney=100;
	FluxSol::FeGrid<2> grid(1.,1.,1.,nex,ney,1);
	
	// FluxSol::FeGrid<2> grid;
	// grid.Create_test(1.,1.,1.,2,2,1);
	// cout << "Total Grid Nodes" <<grid.NumNodes()<<endl;
			// for (int n = 0; n < grid.NumNodes(); n++)
				// cout<<grid.Nod(n).Coords()<<endl;
			
	grid.outstr();

	FluxSol::DoFHandler<2> dofhandler(grid);
	cout << "Global DOFs" << dofhandler.NumDoF()<<endl;
	
	dofhandler.outstr();

	cout << "Assemblying matrix" <<endl;
	
	FluxSol::GaussFullMatrices J,B,dHdrs;
	//FluxSol::ShapeFunctionGroup shfngr = e.CreateShapeFunctionGroup();

	FluxSol::Matrix<double> Kel(8, 8);
	FluxSol::Matrix<double> c(3, 3);

	PETSC_Solver<double,2> solver(dofhandler.NumDoF());

	double E = 200.0e9;
	double nu = 0.33;

	//Plain Stress
	double ck = E / (1 - nu*nu);

	ck = E*(1.-nu)/((1.+nu)*(1.-2*nu));
	c[0][0] = c[1][1] = ck;
	c[0][1] = c[1][0] = ck*nu/(1.-nu);
	c[2][2] = ck*(1 - 2*nu) / (2.*(1.-nu));
	

	for (int e=0;e<grid.NumElem();e++)
	{
		//cout <<"Element: " << e <<endl;
		// TO BE MODIFIED: DONT INCLUDE GRID
        vector <unsigned int> vn = dofhandler.get_dof_indices(e);
		FluxSol::FEValues<2> fev(grid.Elem(e),grid);
		
        J = fev.Jacobian();
		B = fev.shape_grad_matrix();
		
		//B.outstr();
        Kel=0.;
		//cout << "Creating Elemental Stiffness Matrix, GaussPoints: "<< intsch.NumPoints() <<endl;
        for (int g = 0; g < intsch.NumPoints(); g++)
        {
			// cout << "J matrix " << J.Mat(g).outstr()<<endl;
			// cout << "B matrix " << B.Mat(g).outstr()<<endl;
			
            FluxSol::Matrix<double> Kt = B.Mat(g).Tr()*c*B.Mat(g);
            for (int r = 0; r < 8; r++)
				for (int c = 0; c < 8; c++)
					Kel[r][c] += Kt[r][c]*intsch[g].w()*J.Mat(g).det();
			
			//cout <<Kt.outstr();
        }
		

				
        for (int row = 0; row < 8; row++){
            for (int col = 0; col < 8; col++){			
				solver.AddMatVal(vn[row],vn[col],Kel[row][col]);
            }

        }

	}
	
	solver.Flush();
	
	cout << "Applying BCs..."<<endl;
	
	int dof;
	int fix[]={0,1,2*(nex+1)-1};
	for (int f=0;f<3;f++){
		dof=fix[f];
		//cout << "adj dofs"<<endl;
		for (int i = 0; i < dofhandler.Adj_DoF_Number()[dof]; i++){
			int adjdof = dofhandler.AdjDoF(dof,i);
			cout << adjdof << " ";
			//Have to set up only the nonzero matrix values
			solver.SetMatVal(dof,adjdof,0.);
			solver.SetMatVal(adjdof,dof,0.);}
		solver.SetMatVal(dof,dof,1.);	}	

	//2*(nex+1)*ny			
	solver.SetbValues(2*(nex+1)*ney,1000.0);
	

	solver.ViewInfo();
	solver.Solve();
	solver.ViewInfo();
	
	logfile.close();
	
	std::vector<Vec3D> res;
	Vec3D v;
	for (int n = 0; n < grid.NumNodes(); n++){
		for (int d = 0; d < 2; d++)
		    v[d]=solver.X(2*n+d);
		res.push_back(v);
		//cout << v <<endl;
	}
	
	outfile of("output.vtk",grid,res);
	
	logfile << "Stiffness Matrix \n\n";
	logfile << Kel.outstr();

	
	logfile.close();
	return 0;
}

//VTK XML TYPE
//There is also an ASCII type
outfile::outfile(string name, FeGrid<2> &g,const vector<Vec3D> &vres)
:grid(g)
{
	int Rank=0;
	//grid=g;
	string fileName=name;
	ofstream file;
	file.open((fileName).c_str(),ios::out);
	file << "<?xml version=\"1.0\"?>" << endl;
	file << "<VTKFile type=\"UnstructuredGrid\">" << endl;
	file << "<UnstructuredGrid>" << endl;
	file << "<Piece NumberOfPoints=\"" << grid.NumNodes() << "\" NumberOfCells=\"" << grid.NumElem() << "\">" << endl;
	file << "<Points>" << endl;
	file << "<DataArray NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\" >" << endl;
	for (int n=0;n<grid.NumNodes();++n) {
		for (int i=0; i<3; ++i) file<< setw(16) << setprecision(8) << scientific << grid.Nod(n).Coords()[i] << endl;
	}
	file << "</DataArray>" << endl;
	file << "</Points>" << endl;
	file << "<Cells>" << endl;
	
	file << "<DataArray Name=\"connectivity\" type=\"Int32\" format=\"ascii\" >" << endl;
	for (int c=0;c<grid.NumElem();++c) {
		for (int n=0;n<grid.Elem(c).NumNodes();++n) {
			file << grid.Elem(c).NodePos(n) << "\t";}
		file << endl;}
	
	file << "</DataArray>" << endl;
	
	file << "<DataArray Name=\"offsets\" type=\"Int32\" format=\"ascii\" >" << endl;
	int offset=0;
	for (int c=0;c<grid.NumElem();++c) {
		offset+=grid.Elem(c).NumNodes();
		file << offset << endl;}
	file << "</DataArray>" << endl;
	
	file << "<DataArray Name=\"types\" type=\"UInt8\" format=\"ascii\" >" << endl;
	for (int c=0;c<grid.NumElem();++c) {
		// if (grid.Cell(c).Num_Vertex()==4) 
		file << "10" << endl; // Tetra
		// if (grid.Cell(c).Num_Vertex()==8) file << "12" << endl; // Hexa
		// if (grid.Cell(c).Num_Vertex()==6) file << "13" << endl; // Prism
		// if (grid.Cell(c).Num_Vertex()==5) file << "14" << endl; // Pyramid (Wedge)
	}
	file << endl;
	file << "</DataArray>" << endl;;
	
	file << "</Cells>" << endl;
	
	file << "<PointData Scalars=\"scalars\" format=\"ascii\">" << endl;
	
	//Begin data field output
	file << "<DataArray Name=\"";

	file << "x-displacement";
	file << "\" type=\"Float32\" format=\"ascii\" >" << endl;

	//If scalars
	for (int n=0;n<grid.NumNodes();n++)
		    file << vres[n][0] << endl;

	file << "</DataArray>" << endl;

	file << "<DataArray Name=\"";

	file << "y-displacement";
	file << "\" type=\"Float32\" format=\"ascii\" >" << endl;

	//If scalars
	for (int n=0;n<grid.NumNodes();n++)
		    file << vres[n][1] << endl;

	file << "</DataArray>" << endl;

	// // End of data output
	 file << "</PointData>" << endl;
	
	file << "</Piece>" << endl;
	file << "</UnstructuredGrid>" << endl;
	file << "</VTKFile>" << endl;
	file.close();

}