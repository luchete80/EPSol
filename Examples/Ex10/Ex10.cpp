#include "EPSol.h"

//Sparse libraries


using namespace FluxSol;
using namespace std;


int main(int argc,char **args)
{
	cout << "argc " << argc<<endl;
	if (argc>1)
	{

		string filename = args[1];

		cout <<"Opening file..."<< filename <<endl;
		SingleInputFile input(filename);

        //Echoes input file
		//input.ShowData();

		ofstream logfile;
		logfile.open("Logfile.txt");


		const FluxSol::FEIntegrationScheme intsch;
		FluxSol::GaussMatrices H;

		FluxSol::GaussFullMatrices J,B,dHdrs;
		FluxSol::ShapeFunctionGroup shfngr;



		//Plain Strain
		//ck = E*(1.-nu)/((1.+nu)*(1.-2*nu));
		//c[0][0] = c[1][1] = ck;
		//c[0][1] = c[1][0] = ck*nu/(1.-nu);
		//c[2][2] = ck*(1 - 2*nu) / (2.*(1.-nu));

		//C = E / (1 - nu*nu)*[1 nu 0; nu 1 0; 0 0 (1 - nu) / 2];

		int totrows;
		int row,col;

		// INIT MESH

		const int dim=2;

		//Mesh
		FluxSol::FeGrid<dim> grid(input.Grid());
		cout << "Total Grid Nodes" <<grid.NumNodes()<<endl;
		FluxSol::DoFHandler<dim> dofhandler(grid);
		//cout << dofhandler.outstr();

		//MESH INFO
		//cout << "Mesh Info" <<endl;
		//for (int e=0;e<grid.NumElem();e++)
		//	for (int n=0;n<4;n++)
		//		cout<<grid.Elem(e).NodePos(n)<<endl;


		//Dummy assembly Matrix
		//cout << "Creating Matrix... " << dofhandler.NumDoF()<<endl;
		//FluxSol::Matrix<double> Kgi(dofhandler.NumDoF(),dofhandler.NumDoF());

		PETSC_Solver<double,dim> solver(dofhandler.NumDoF());
		//TO MODIFY
		//Can be an integer or a non constant vector
		//solver.PreAllocateRows(dofhandler.Adj_DoF_Number());

		//cout << "test" <<endl;
		//Q4PS elpru;
		//FluxSol::FEValues<dim> ff(elpru);

		//c=ff.MatMatrix();

		//Assemble of Matrix

		//TO MODIFY: ASSUMING THAT FIELD_DIM IS CONSTANT
		const int field_dim=grid.Elem(0).Field_Dim();

        FluxSol::Matrix<double> Kel(grid.Elem(0).DOF(), grid.Elem(0).DOF());
		FluxSol::Matrix<double> c;


		cout << "Assemblying matrix" <<endl;
		for (int e=0;e<grid.NumElem();e++)
		{
			//cout <<"Element: " << e <<endl;
			// TO BE MODIFIED: DONT INCLUDE GRID
			//FluxSol::FEValues<2> fev(grid.Elem(e), grid);
			vector <unsigned int> vn = dofhandler.get_dof_indices(e);
			cout << "dof index initiated" <<endl;
			//TO MODIFY
			grid.Elem(e).AssignMat(input.Mat());	//TO MODIFY, MATERIAL VECTOR
			//for (int i=0;i<vn.size();i++)
			//	cout<<"GLOBAL DOFS"<<vn[i]<<endl;

			cout << "Material initialized ..."<<endl;
			//FluxSol::FEValues<2> fev(qel,g);
			FluxSol::FEValues<dim> fev(grid.Elem(e),grid);
			cout << "Finite Element Values initiated ..."<<endl;
			J = fev.Jacobian();
			B = fev.shape_grad_matrix();

			cout << "Element Field dim " << grid.Elem(e).Field_Dim() <<endl;

			cout << "Initializing Material Matrix..."<<endl;
			c=fev.MatMatrix();

			cout << "Matrix initialized..." <<endl;
			cout << "ELEMENT MATRIX " << c.outstr() <<endl;

			Kel.Clear();
			cout << "Creating Elemental Stiffness Matrix, GaussPoints: "<< intsch.NumPoints() <<endl;
			for (int g = 0; g < intsch.NumPoints(); g++)
			{

				FluxSol::Matrix<double> Kt = B.Mat(g).Tr()*c*B.Mat(g);

				//Matrix <double> Kt =
				//Look throug Element DOFS
				for (int r = 0; r < grid.Elem(e).DOF(); r++)
					for (int c = 0; c < grid.Elem(e).DOF(); c++)
						Kel[r][c] += Kt[r][c]*intsch[g].w()*J.Mat(g).det();

				//cout <<Kt.outstr();
			}

			//logfile << "Stiffness Matrix \n\n";
			//cout << Kel.outstr();

			//cout <<"Dof Indices"<<endl;
			//for(int row = 0; row < grid.Elem(e).DOF(); row++)
			//	cout <<vn[row]<<endl;

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
					//Kgi[vn[row]][vn[col]]+=Kel[row][col];

					solver.AddMatVal(vn[row],vn[col],Kel[row][col]);
				}

			}

		}

		//cout<<Kgi.outstr();

		//////////////////////////////////
		//Applying Boundary Condition ////
		//////////////////////////////////
		solver.Flush();
		int node,dir,dof;
		cout << "BFix Number: " <<input.BFixNumber() << endl;

		for (int bcdof=0;bcdof<input.BFixNumber();bcdof++)
		{
			cout<<input.BFix(bcdof).Show();

			node=input.BFix(bcdof).NodeId();
			std::vector<double> val=input.BFix(bcdof).Value();

			const vector<bool>fixed=input.BFix(bcdof).Fix();
			//TO MODIFY, BETTER WITH A SINGLE VECTOR CONTAINING THE ONLY FIZED DIR
			for (int d=0;d<field_dim;d++)
			{
				cout << "dimension " <<d<<endl;
				if (fixed[d])
				{
					//cout << "node " <<node << " dim " << d << " fixed"<<endl;
					//If is a original nonzero value, must be added
					dof=field_dim*node+d;
                    solver.ApplyBCOnDoF(dof, dofhandler, val[d]);
				}
			}

		}

		//Fixing y=0,
		cout << "Applying fixities..."<<endl;

		for (int ii=0;ii<grid.NumNodes();ii++)
        {
            Vec3D v=grid.Nod(ii).Coords();
            if (v.y()==0.0)
            {
                int dof=field_dim*dofhandler.Vertex_DoF(ii);
                solver.ApplyBCOnDoF(dof, dofhandler, 0.0);
            }
        }

        ///////////////////////////////////
        //APPLYING CONVECTION BOUNDARIES //
        double hc=200.0; //W/m2-K
        double To=20.0;
        double rval=-hc*To;
        cout << "Applying Convection..."<<endl;
        for (int ii=0;ii<grid.NumNodes();ii++)
        {
            int dof=field_dim*dofhandler.Vertex_DoF(ii);
            solver.SetbValues(dof, rval);
            solver.AddMatVal(dof,dof,-hc);
        }

		cout << "BLoad Number: " <<input.BLoadNumber() << endl;


        //DOF Order Against Node Order
        //Is directly assigned??
		for (int bcdof=0;bcdof<input.BLoadNumber();bcdof++)
		{
			//cout<<input.BLoad(bcdof).Show();

			node=input.BLoad(bcdof).NodeId();

			const vector<double>load=input.BLoad(bcdof).Load();
			//TO MODIFY, BETTER WITH A SINGLE VECTOR CONTAINING THE ONLY FIZED DIR
			for (int d=0;d<field_dim;d++)
			{
					dof=field_dim*node+d;
                    cout << "Applying dof "<<dof<<endl;
					solver.SetbValues(dof,load[d]);

					//cout << "DOF / LOAD" << dof << " " <<load[d] <<endl;
			}

		}

		cout << "Applying loads ..."<<endl;
		for (int ii=0;ii<grid.NumNodes();ii++)
        {
            Vec3D v=grid.Nod(ii).Coords();
            if (v.y()==0.004 && v.x()<=0.006)
            {
                cout << "Node " <<ii<<endl;
                int dof=field_dim*dofhandler.Vertex_DoF(ii);
                double load;
                solver.SetbValues(dof, load);
            }
        }


		solver.ViewInfo();
		solver.Solve();
		solver.ViewInfo();

        //Color nodes
		//Useful Nodes

        ofstream p3dfile;
		p3dfile.open("p3dfile.txt");
        p3dfile << "x y z temp"<<endl;
		for (int ii=0;ii<grid.NumNodes();ii++)
        {
            Vec3D v=grid.Nod(ii).Coords();
            p3dfile <<v.x()<< " " << v.y() << " " << v.z() << " ";
            for (int d=0;d<field_dim;d++)
            {
                int dof=field_dim*dofhandler.Vertex_DoF(ii)+d;
                p3dfile <<solver.X(dof) << " " ;
            }
            p3dfile << endl;
        }

	}//Enf if arguments
	else
	{
		cout << "No input file. See epsol command syntax."<<endl;
	}

	return 0;
}
