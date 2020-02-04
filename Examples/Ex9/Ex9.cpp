#include "EPSol.h"

//Sparse libraries


using namespace FluxSol;
using namespace std;

///////////////////////////
//TRANSIENT THERMAL EXAMPLE
///////////////////////////

int main(int argc,char **args)
{
	cout << "argc " << argc<<endl;
	if (argc>1)
	{

		string filename = args[1];

		cout <<"Opening file..."<< filename <<endl;
		SingleInputFile input(filename);

		input.ShowData();

		ofstream logfile;
		logfile.open("Logfile.txt");


		const FluxSol::FEIntegrationScheme intsch;
		FluxSol::GaussFullMatrices H;

		FluxSol::GaussFullMatrices J,B,dHdrs;
		FluxSol::ShapeFunctionGroup shfngr;


		FluxSol::Matrix<double> Kel(8, 8);
		FluxSol::Matrix<double> Mel(8, 8);
		FluxSol::Matrix<double> c(3, 3);


		int totrows;
		int row,col;
		int tottime=10.;

		double rho, cp;

		rho=cp=1.0;

		//INTEGRATION DATA
		double dt=1.0;
		double theta=0.5;

		// INIT MESH

		const int dim=2;

		//Mesh
		FluxSol::FeGrid<dim> grid(input.Grid());
		cout << "Total Grid Nodes" <<grid.NumNodes()<<endl;
		FluxSol::DoFHandler<dim> dofhandler(grid);
		cout << dofhandler.outstr();

		//MESH INFO
		cout << "Mesh Info" <<endl;
		for (int e=0;e<grid.NumElem();e++)
			for (int n=0;n<4;n++)
				cout<<grid.Elem(e).NodePos(n)<<endl;

		//Dummy assembly Matrix
		FluxSol::Matrix<double> Kgi(dofhandler.NumDoF(),dofhandler.NumDoF());

		PETSC_Solver<double,dim> solver(dofhandler.NumDoF());
		//TO MODIFY
		//Can be an integer or a non constant vector
		//solver.PreAllocateRows(dofhandler.Adj_DoF_Number());

		//Force elemental Vector History, previous time
		cout << "Force Vectors" <<endl;
		Vector<double> Ft_el(grid.Elem(0).DOF());
		Vector<double> Tt_el(grid.Elem(0).DOF());

		Vector<double> F_el(grid.Elem(0).DOF());

		//////////////////////
		// TIME INTEGRATION //
		//////////////////////
		// EQN
		// (C/dt + K * theta) T_(t+dT) = (1-theta) * F_t + theta * F_(t+dT) + (C/dT -K*(1-theta) ) * T_t
		cout << "Time Loop ..."<<endl;
		for (double t=0.0;t<tottime;t+=dt)
		{
			//TO MODIFY: ASSUMING THAT FIELD_DIM IS CONSTANT
			const int field_dim=grid.Elem(0).Field_Dim();

			cout << "Assemblying matrix" <<endl;
			for (int e=0;e<grid.NumElem();e++)
			{
				// --------------------- GET FORCE VECTORS
				//cout << "Getting force vectors"<<endl;
				//SAVE History for time integration
				//TO ANALYZE: COULD BE GET VIA Solver::B() function the entire force assembled vector
				//IT IS CONVENIENT
				vector <unsigned int> vn = dofhandler.get_dof_indices(e);
				for (int idof=0;idof<grid.Elem(e).DOF();idof++)
				{
					int gdof=vn[idof]; //global dof
					cout << "Value" << solver.B(gdof);
					Ft_el[idof]=solver.B(gdof);
					Tt_el[idof]=solver.X(gdof);
				}


				//cout <<"Element: " << e <<endl;
				// TO BE MODIFIED: DONT INCLUDE GRID
				//FluxSol::FEValues<2> fev(grid.Elem(e), grid);
				vector <unsigned int> vn = dofhandler.get_dof_indices(e);
				//cout << "dof index initiated" <<endl;
				//TO MODIFY
				grid.Elem(e).AssignMat(input.Mat());	//TO MODIFY, MATERIAL VECTOR
				//for (int i=0;i<vn.size();i++)
				//	cout<<"GLOBAL DOFS"<<vn[i]<<endl;

				//cout << "Material initialized ..."<<endl;
				//FluxSol::FEValues<2> fev(qel,g);
				FluxSol::FEValues<dim> fev(grid.Elem(e),grid);
				//cout << "Finite Element Values initiated ..."<<endl;
				J = fev.Jacobian();
				H = fev.shape_value_matrix();
				B = fev.shape_grad_matrix();

				//cout << "Element Field dim " << grid.Elem(e).Field_Dim() <<endl;

				cout << "Initializing Material Matrix..."<<endl;
				c=fev.MatMatrix();

				//cout << "Matrix initialized..." <<endl;
				//cout << "ELEMENT MATRIX " << c.outstr() <<endl;

				Kel.Clear();
				Mel.Clear();

				//cout << "Creating Elemental Stiffness Matrix, GaussPoints: "<< intsch.NumPoints() <<endl;

				for (int g = 0; g < intsch.NumPoints(); g++)
				{

					FluxSol::Matrix<double> Kt = B.Mat(g).Tr()*c*B.Mat(g);

					FluxSol::Matrix<double> Mm = H.Mat(g).Tr()*H.Mat(g)*rho*cp;		//Mass Matrix

					//Matrix <double> Kt =
					//Look throug Element DOFS
					//TO MODIFY, SIMPLIFY INTEGRATION TRHOUGH GAUSS POINTS
					for (int r = 0; r < grid.Elem(e).DOF(); r++)
						for (int c = 0; c < grid.Elem(e).DOF(); c++)
						{
							Kel[r][c] += Kt[r][c]*intsch[g].w()*J.Mat(g).det();
							Mel[r][c] += Mm[r][c]*intsch[g].w()*J.Mat(g).det();
						}

					//cout <<Kt.outstr();
				}//Gauss Points

				//logfile << "Stiffness Matrix \n\n";
				//logfile << Kel.outstr();

				//cout <<"Dof Indices"<<endl;
				//for(int row = 0; row < grid.Elem(e).DOF(); row++)
				//	cout <<vn[row]<<endl;


				//FORCE VECTOR
				F_el=Ft_el*(1.-theta)+(Mel*(1./dt)-Kel*theta)*Tt_el;

				for (int row = 0; row < grid.Elem(e).DOF(); row++)
				{
					for (int col = 0; col < grid.Elem(e).DOF(); col++)
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

						// (C/dt + K * theta) T_(t+dT) = (1-theta) * F_t + theta * F_(t+dT) + (C/dT -K*(1-theta) ) * T_t
						solver.AddMatVal(vn[row],vn[col],Mel[row][col]/dt+Kel[row][col]*theta);
					}
					//..= (1-theta) * F_t + theta * F_(t+dT) + (C/dT -K*(1-theta) ) * T_t
					solver.SetbValues(vn[row],(1.-theta)*);
				}

			}//Loop Through Elements
			///////////////////////////////////////////////////////////////////////////////////////


			//cout<<Kgi.outstr();

			//////////////////////////////////
			//Applying Boundary Condition ////
			//////////////////////////////////
			solver.Flush();
			int node,dir,dof;
			//cout << "BFix Number: " <<input.BFixNumber() << endl;

			for (int bcdof=0;bcdof<input.BFixNumber();bcdof++)
			{
				cout<<input.BFix(bcdof).Show();

				node=input.BFix(bcdof).NodeId();

				const vector<bool>fixed=input.BFix(bcdof).Fix();
				//TO MODIFY, BETTER WITH A SINGLE VECTOR CONTAINING THE ONLY FIZED DIR
				for (int d=0;d<field_dim;d++)
				{
					cout << "dimension" <<endl;
					if (fixed[d])
					{
						cout << "node " <<node << " dim " << d << " fixed"<<endl;
						//If is a original nonzero value, must be added
						dof=field_dim*node+d;

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

			cout << "BLoad Number: " <<input.BLoadNumber() << endl;

			for (int bcdof=0;bcdof<input.BLoadNumber();bcdof++)
			{
				//cout<<input.BLoad(bcdof).Show();

				node=input.BLoad(bcdof).NodeId();

				const vector<double>load=input.BLoad(bcdof).Load();
				//TO MODIFY, BETTER WITH A SINGLE VECTOR CONTAINING THE ONLY FIZED DIR
				for (int d=0;d<field_dim;d++)
				{
						dof=field_dim*node+d;
						solver.SetbValues(dof,load[d]);

						//cout << "DOF / LOAD" << dof << " " <<load[d] <<endl;
				}

			}


			//solver.ViewInfo();
			solver.Solve();
			//solver.ViewInfo();

		} //TIME LOOP

	}//Enf if arguments
	else
	{
		cout << "No input file. See epsol command syntax."<<endl;
	}

	return 0;
}
