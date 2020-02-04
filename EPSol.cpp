#include "EPSol.h"

//Sparse libraries


using namespace FluxSol;
using namespace std;

void Ex11();
int Ex9(const string &cad);
int Ex12(const string &cad);

int main(int argc,char **args)
{
	cout << "argc " << argc<<endl;
	//if (argc>1)
	//if (0)
	//if (argc>1)
	if (0)
	{
		const int dim=3;



		//cout <<"Opening file..."<< filename <<endl;
		//string filename = args[1];
		//SingleInputFile<dim> input(filename);

		SingleInputFile<dim> input("1Elem-TH.in");

        //Echoes input file
		//input.ShowData();

		ofstream logfile;
		logfile.open("Logfile.txt");


		FluxSol::GaussFullMatrices J,B,H,dHdrs;
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




        //double val = shfngr.ShapeFn(0).Val(0.577, 0.577, 0.);
        shfngr = grid.Elem(e).CreateShapeFunctionGroup();
        logfile << shfngr.ShapeFn(2).outstr();

        for (int d = 0; d < 8; d++)
        {
            cout << "\n Shape Fn " << d << "\n";
            vector<FluxSol::ShapeFunction> df = shfngr.ShapeFn(d).diff();
            ShapeFunction f = shfngr.ShapeFn(d);
            cout << "F valuated"<<endl;
            cout << f.outstr()<<endl;

            cout << "\n Local Diff r Coeff \n";
            cout << df[0].outstr();
            cout << "\n Local Diff s Coeff \n";
            cout << df[1].outstr();
            cout << "\n Local Diff t Coeff \n";
            cout << df[2].outstr();
        }





			cout << "Material initialized ..."<<endl;
			//FluxSol::FEValues<2> fev(qel,g);
			FluxSol::FEValues<dim> fev(grid.Elem(e),grid);
			cout << "Finite Element Values initiated ..."<<endl;
			J = fev.Jacobian();
			B = fev.shape_grad_matrix();
			H = fev.shape_value_matrix();


			cout << "Element Field dim " << grid.Elem(e).Field_Dim() <<endl;

			cout << "Initializing Material Matrix..."<<endl;
			c=fev.MatMatrix();

			cout << "Matrix initialized..." <<endl;
			cout << "ELEMENT MATRIX " << c.outstr() <<endl;

			Kel.Clear();
			FEIntegrationScheme intsch(grid.Elem(e).GaussOrder(),dim);

			cout << "Creating Elemental Stiffness Matrix, GaussPoints: "<< intsch.NumPoints() <<endl;
			for (int g = 0; g < intsch.NumPoints(); g++)
			{

				FluxSol::Matrix<double> Kt = B.Mat(g).Tr()*c*B.Mat(g);

				//Matrix <double> Kt =
				//Look throug Element DOFS
				for (int r = 0; r < grid.Elem(e).DOF(); r++)
					for (int c = 0; c < grid.Elem(e).DOF(); c++)
						Kel[r][c] += Kt[r][c]*intsch[g].w()*J.Mat(g).det();
                cout << "det" << J.Mat(g).det()<<endl;
                cout << "det_" << J.Mat(g).det_()<<endl;
                cout << "J" << J.Mat(g).outstr()<<endl;
                //cout << "Gauss Point Weight "<< intsch[g].w() <<endl;
                //cout << "Gauss Coords" << intsch[g].Coords().Imprimir_Coord()<<endl;

                cout << "H valuated " <<H.Mat(g).outstr();
                cout << "B valuated " <<B.Mat(g).outstr();
				//cout <<Kt.outstr();
			}

			logfile << "Stiffness Matrix \n\n";
			cout << Kel.outstr();

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
        //Applying Global Boundary Condition ////
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


		cout << "BLoad Number: " <<input.BLoadNumber() << endl;


        //DOF Order Against Node Order
        //Is directly assigned??
        solver.SetbValues(0.);
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

					cout << "DOF / LOAD" << dof << " " <<load[d] <<endl;
			}

		}

		solver.SetbValues(7,1.0);


		//solver.ViewInfo();
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


    //Ex9(argc,args);

    Ex9("1Elem-TH.in"); //TRANSIENT
    //Ex11();
	return 0;
}

///////////////////
//Transient Example
int Ex9(const string &cad)
{
        const int dim=3;

		string filename = cad;

		cout <<"Opening file..."<< filename <<endl;
		SingleInputFile<dim> input(filename);

		input.ShowData();

		ofstream logfile;
		logfile.open("Logfile.txt");



		//const FluxSol::FEIntegrationScheme intsch;
		FluxSol::GaussFullMatrices H;

		FluxSol::GaussFullMatrices J,B,dHdrs;
		FluxSol::ShapeFunctionGroup shfngr;


		FluxSol::Matrix<double> c(3, 3);


		int totrows;
		int row,col;
		int tottime=50.;

		double rho, cp;

		rho=cp=1.0;

		//INTEGRATION DATA
		double dt=1.0;
		double theta=1.0;

		// INIT MESH


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
		Vector<double> Ftpdt_el(grid.Elem(0).DOF());

		Vector<double> F_el(grid.Elem(0).DOF());
		Vector<double> R_el(grid.Elem(0).DOF());


        //Load vectors
        //std::vector<Vector <double> > Fel_t;
        //Fel_t.assign(grid.NumElem(),grid.Elem(0).DOF());

        Vector<double> F_glob_tind(dofhandler.NumDoF()); //Global load vector, time independent

		FluxSol::Matrix<double> Kel(grid.Elem(0).DOF(), grid.Elem(0).DOF());
		FluxSol::Matrix<double> Mel(grid.Elem(0).DOF(), grid.Elem(0).DOF());

        const int field_dim=grid.Elem(0).Field_Dim();
		// TIME INDEPENDENT //
        //////////////////////////
        //APPLYING GLOBAL LOADS //
        //////////////////////////
        cout << "BLoad Number: " <<input.BLoadNumber() << endl;

        for (int bcdof=0;bcdof<input.BLoadNumber();bcdof++)
        {
            //cout<<input.BLoad(bcdof).Show();

            int node=input.BLoad(bcdof).NodeId();

            const vector<double>load=input.BLoad(bcdof).Load();
            //TO MODIFY, BETTER WITH A SINGLE VECTOR CONTAINING THE ONLY FIZED DIR
            for (int d=0;d<field_dim;d++)
            {
                    int dof=field_dim*node+d;
                    //solver.AddbValues(dof,theta*load[d]);
                    F_glob_tind[dof]=load[d];
                    //cout << "DOF / LOAD" << dof << " " <<load[d] <<endl;
            }

        }

		//////////////////////
		// TIME INTEGRATION //
		//////////////////////
		// EQN
		// (C/dt + K * theta) T_(t+dT) = (1-theta) * F_t + theta * F_(t+dT) + (C/dT -K*(1-theta) ) * T_t
		cout << "Time Loop ..."<<endl;
		for (double t=0.0;t<tottime;t+=dt)
		{
			//TO MODIFY: ASSUMING THAT FIELD_DIM IS CONSTANT


            //Save elemental loads


            //////////////////////////////////////////////////////
            ///////// APPLYING ELEMENTAL NODES (IF ANY) //////////


			////////////////////////////////////
			///// LOOP THROUGH ELEMENTS ////////
			////////////////////////////////////


			cout << "Assemblying matrix" <<endl;
			for (int e=0;e<grid.NumElem();e++)
			{
			    FEIntegrationScheme intsch(grid.Elem(e).GaussOrder(),dim);
				// --------------------- GET FORCE VECTORS
				//cout << "Getting force vectors"<<endl;
				//SAVE History for time integration
				//TO ANALYZE: COULD BE GET VIA Solver::B() function the entire force assembled vector
				//IT IS CONVENIENT
				vector <unsigned int> vn = dofhandler.get_dof_indices(e);
				for (int idof=0;idof<grid.Elem(e).DOF();idof++)
				{
					int gdof=vn[idof]; //global dof
					Tt_el[idof]=solver.X(gdof); //Solution Loads Previous time
					Ftpdt_el[idof]=F_glob_tind[gdof]; //Solution Loads Previous time
				}

				if (t==0.)   Ft_el=0.;
                else
                {
                    for (int idof=0;idof<grid.Elem(e).DOF();idof++)
                    {
                        int gdof=vn[idof]; //global dof
                        Ft_el[idof]=F_glob_tind[gdof]; //Solution Loads Previous time
                    }
                }

				cout << "Element Loads"<<endl;
				cout << Ft_el.outstr()<<endl;

				cout << "Element Temperatures"<<endl;
				cout << Tt_el.outstr()<<endl;

				//cout <<"Element: " << e <<endl;
				// TO BE MODIFIED: DONT INCLUDE GRID
				//FluxSol::FEValues<2> fev(grid.Elem(e), grid);
				vector <unsigned int> vne = dofhandler.get_dof_indices(e);
				cout <<"dof assigned"<<endl;
				//cout << "dof index initiated" <<endl;
				//TO MODIFY
				grid.Elem(e).AssignMat(input.Mat());	//TO MODIFY, MATERIAL VECTOR
				//for (int i=0;i<vn.size();i++)
				//	cout<<"GLOBAL DOFS"<<vn[i]<<endl;

				//cout << "Material initialized ..."<<endl;
				//FluxSol::FEValues<2> fev(qel,g);
				cout << "Creating FEValues"<<endl;
				FluxSol::FEValues<dim> fev(grid.Elem(e),grid);
				cout << "Fev Created"<<endl;
				//cout << "Finite Element Values initiated ..."<<endl;
				J = fev.Jacobian();
				H = fev.shape_value_matrix();
				cout << "H created"<<endl;
				B = fev.shape_grad_matrix();

				//cout << "Element Field dim " << grid.Elem(e).Field_Dim() <<endl;

				cout << "Initializing Material Matrix..."<<endl;
				c=fev.MatMatrix();

				//cout << "Matrix initialized..." <<endl;
				cout << "ELEMENT Material MATRIX " << c.outstr() <<endl;

				Kel.Clear();
				Mel.Clear();

				cout << "Creating Elemental Stiffness Matrix, GaussPoints: "<< intsch.NumPoints() <<endl;

				for (int g = 0; g < intsch.NumPoints(); g++)
				{

                    //cout << "gauss point " << g<<endl;
					FluxSol::Matrix<double> Kt = B.Mat(g).Tr()*c*B.Mat(g);

					FluxSol::Matrix<double> Mm = H.Mat(g).Tr()*H.Mat(g)*rho*cp;		//Mass Matrix
					cout << "H valuated " <<H.Mat(g).outstr();
					cout << "B valuated " <<B.Mat(g).outstr();
                    cout << "Determinant" <<J.Mat(g).det();
                    cout << "Gauss Point" <<intsch[g].Coords().x()<< " " <<
                                            intsch[g].Coords().y() << " " <<
                                            intsch[g].Coords().z() << " " << endl;
					//Matrix <double> Kt =
					//Look throug Element DOFS
					//TO MODIFY, SIMPLIFY INTEGRATION TRHOUGH GAUSS POINTS

					for (int r = 0; r < grid.Elem(e).DOF(); r++)
						for (int c = 0; c < grid.Elem(e).DOF(); c++)
						{
						    //cout << "r c "<<r<< " " << c<<endl;
							Kel[r][c] += Kt[r][c]*intsch[g].w()*J.Mat(g).det();
							Mel[r][c] += Mm[r][c]*intsch[g].w()*J.Mat(g).det();
						}



					//cout <<Kt.outstr();
				}//Gauss Points

				cout << "Stiffness Matrix \n\n";
				cout << Kel.outstr();

				cout << "Mass Matrix \n\n";
				cout << Mel.outstr();


				//cout <<"Dof Indices"<<endl;
				//for(int row = 0; row < grid.Elem(e).DOF(); row++)
				//	cout <<vn[row]<<endl;


                //Get previous Elemental Load Force vector



                //Get elemental force vector



				//ELEMENTAL FORCE VECTOR
				//..= (1-theta) * F_t + theta * F_(t+dT) + (C/dT -K*(1-theta) ) * T_t
				//To this expression misses out the t+dt term affected by alpha, which is nodal
				Vector <double> temp=(Ft_el*0.5);
				cout << "0.5 Ft_el " <<temp.outstr()<<endl;
				F_el=(Ft_el*(1.0-theta))+(Ftpdt_el*theta)+((Mel*(1./dt)-Kel*(1.0-theta))*Tt_el);
				cout << "R el" << F_el.outstr()<<endl;
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
						Kgi[vne[row]][vn[col]]+=Kel[row][col];

						// (C/dt + K * theta) T_(t+dT) = (1-theta) * F_t + theta * F_(t+dT) + (C/dT -K*(1-theta) ) * T_t
						solver.AddMatVal(vne[row],vne[col],Mel[row][col]/dt+Kel[row][col]*theta);
					}
					//HERE ADD THE ELEMENTAL CONTRIBUTION DUE TO TRANSIENT
					//..= (1-theta) * F_t + theta * F_(t+dT) + (C/dT -K*(1-theta) ) * T_t
					solver.AddbValues(vne[row],F_el[row]);
				}

			}//Loop Through Elements
			///////////////////////////////////////////////////////////////////////////////////////


			//cout<<Kgi.outstr();

			//////////////////////////////////
			//Applying Global Boundary Condition ////
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



			//solver.ViewInfo();
			solver.Solve();
			solver.ViewInfo();

            //RESETTING GLOBAL MATRICES
            //Clear Vector and Matrix Solver Values
            solver.ResetMatrix(dofhandler);
            solver.ResetB(dofhandler);      //Global Force Vector, A.K.A. R

		} //TIME LOOP



	return 0;
}

//////////////////////////
//Surface Loads Example //

void Ex11()
{
	ofstream logfile;
	logfile.open("Logfile.txt");

	//Element can be construted from vertex too
	//std::_Vertex
    cout << "Defining Nodes"<<endl;
	std::vector<Node> v;
    cout << "pb"<<endl;
	//TO MODIFY: "ADDNODEFUNCTION IN GRID"
	v.push_back(Node(0, 1.0, 1.0, 0.0));
	cout << "pb1"<<endl;
	v.push_back(Node(1, 0.0, 1.0, 0.0));
	v.push_back(Node(2, 0.0, 0.0, 0.0));
	v.push_back(Node(3, 1.0, 0.0, 0.0));


	//const FEIntegrationScheme intsch;

    cout << "Creating element..."<<endl;
	Q4TH<2> e(v);
	e.Set_Nodes(4, 0, 1, 2, 3);
	logfile << "Element Gauss Order: " << e.GaussOrder() << "\n\n";

	GaussMatrices H = e.H();

	std::vector< Element<2> > ve;
	ve.push_back(e);
	cout << "Creating Mesh..."<<endl;
	FeGrid <2> g(v, ve);
	FEValues<2> fev(e, g);
}

//3D Element
int Ex12(const string &cad)
{
		string filename = cad;

		cout <<"Opening file..."<< filename <<endl;

        int dim=2;

		if (dim==2)
		SingleInputFile<2> input(filename);
        else if (dim==3)
		SingleInputFile<3> input(filename);

}
