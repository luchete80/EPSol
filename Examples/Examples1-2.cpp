#include "EPSol.h"

//Sparse libraries
#include "laspack.h"

//using namespace FluxSol;
//using namespace std;

//#include <petscksp.h>


int main()
{

	ofstream logfile;
	logfile.open("Logfile.txt");

	//Element can be construted from vertex too
	//std::_Vertex

	std::vector<FluxSol::Node> v;

	//TO MODIFY: "ADDNODEFUNCTION IN GRID"
	v.push_back(FluxSol::Node(0, 1.0, 1.0,0.0));
	v.push_back(FluxSol::Node(1, 0.0, 1.0, 0.0));
	v.push_back(FluxSol::Node(2, 0.0, 0.0, 0.0));
	v.push_back(FluxSol::Node(3, 1.0, 0.0, 0.0));

	const FluxSol::FEIntegrationScheme intsch;


	FluxSol::QuadLinearElement e(v);
	e.Set_Nodes(4,0,1,2,3);
	logfile << "Element Gauss Order: " << e.GaussOrder()<<"\n\n";

	FluxSol::GaussMatrices H = e.H();

	std::vector< FluxSol::Element<2> > ve;
	ve.push_back(e);
	FluxSol::FeGrid <2> g(v, ve);
	FluxSol::FEValues<2> fev(e, g);

	FluxSol::GaussFullMatrices J = fev.Jacobian();

	FluxSol::GaussFullMatrices B = fev.shape_grad_matrix();

	FluxSol::GaussFullMatrices dHdrs = fev.shape_localgrad_matrix();

	FluxSol::ShapeFunctionGroup shfngr = e.CreateShapeFunctionGroup();

	double val= shfngr.ShapeFn(0).Val(0.577, 0.577, 0.);

	logfile << shfngr.ShapeFn(2).outstr();

	for (int d = 0; d < 4; d++)
	{
		logfile << "\n Shape Fn "<<d <<"\n";
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
	c[2][2] = ck*(1-nu)/2.;

	//Plain Strain
	//ck = E*(1.-nu)/((1.+nu)*(1.-2*nu));
	//c[0][0] = c[1][1] = ck;
	//c[0][1] = c[1][0] = ck*nu/(1.-nu);
	//c[2][2] = ck*(1 - 2*nu) / (2.*(1.-nu));

	//C = E / (1 - nu*nu)*[1 nu 0; nu 1 0; 0 0 (1 - nu) / 2];

	for (int g = 0; g < intsch.NumPoints(); g++)
	{
		FluxSol::Matrix<double> Kg = B.Mat(g).Tr()*c*B.Mat(g);
		for (int r = 0; r < 8; r++)
		for (int c = 0; c < 8; c++)
			Kel[r][c] += Kg[r][c]*intsch[g].w()*J.Mat(g).det();

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
			Q_SetEntry(&K, r, c,c+1,Kel[r-1][c]);

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
	double a= U.Cmp[1];

	vector<double> Ui;
	Ui.assign(8,0.);

	logfile << "\n Results \n";
	for (int j = 0; j < totrows; j++)
	{
		Ui[j] = U.Cmp[j + 1];
		logfile << Ui[j] << "\n";
	}


	Q_Destr(&K);

	//Mesh
	FluxSol::FeGrid<2> grid(2.,1.,1.,2,1,1);
	cout << "Total Grid Nodes" <<grid.NumNodes()<<endl;

	FluxSol::DoFHandler<2> dofhandler(grid);


    //Global System
    QMatrix Kg;
	Vector Ug, Rg;

	totrows = dofhandler.NumDoF();
	cout << "Total DOFs" << totrows<<endl;
	Symmetry = False;
	Q_Constr(&Kg, "K", totrows, Symmetry, Rowws, Normal, True);
	V_Constr(&Ug, "U", totrows, Normal, True);
	V_Constr(&Rg, "R", totrows, Normal, True);

	FluxSol::Matrix<double> Kgi(12,12);

    for (int r=1;r<=12;r++)
        Q_SetLen(&Kg, r, 12);

	//Assemble of Matrix
	for (int e=0;e<grid.NumElem();e++)
	{
        vector <unsigned int > vn = dofhandler.get_dof_indices(e);

        for (int i=0;i<8;i++)
            cout<<vn[i]<<endl;

        FluxSol::FEValues<2> feval(grid.Elem(e));

        FluxSol::GaussFullMatrices J = fev.Jacobian();

        FluxSol::GaussFullMatrices B = fev.shape_grad_matrix();

        Kel.Clear();
        for (int g = 0; g < intsch.NumPoints(); g++)
        {
            cout << "Peso" << intsch[g].w();
            FluxSol::Matrix<double> Kt = B.Mat(g).Tr()*c*B.Mat(g);
            for (int r = 0; r < 8; r++)
            for (int c = 0; c < 8; c++)
                Kel[r][c] += Kt[r][c]*intsch[g].w()*J.Mat(g).det();

        }

            logfile << "Stiffness Matrix \n\n";
            logfile << Kel.outstr();

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


        //logfile << "Global Stiffness Matrix \n\n";
        //logfile << K.outstr();

    for (int row = 0; row < 12; row++)
        for (int col = 0; col < 12; col++)
            Q_SetEntry(&Kg, row+1, col,col+1,Kgi[row][col]);


	for (int i = 0; i < 12; i++)
	{
		//Q_SetEntry(&Kg, r from 1, c from 0,c from 1,Kel[r][c]);
		Q_SetEntry(&Kg, 1, i,i+1,0.0);
		Q_SetEntry(&Kg, i+1, 0,1,0.0);

		Q_SetEntry(&Kg, 2, i,i+1,0.0);
		Q_SetEntry(&Kg, i+1, 1,2,0.0);

		Q_SetEntry(&Kg, 6, i,i+1,0.0);
		Q_SetEntry(&Kg, i+1, 5,6,0.0);
	}

	Q_SetEntry(&Kg, 1, 0,1,1.0);
	Q_SetEntry(&Kg, 2, 1,2,1.0);
	Q_SetEntry(&Kg, 6, 5,6,1.0);

    V_SetAllCmp(&Rg, 0.0);
	V_SetCmp(&Rg, 7, 1000.0);


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


	//////////////////// PETSC SOLVER //////////////////////
//	Vec x, b, u; /* approx solution, RHS, exact solution */
//	Mat A; /* linear system matrix */
//	KSP ksp; /* linear solver context */
//	PC pc; /* preconditioner context */
//	PetscReal norm; /* norm of solution error */
//	PetscErrorCode ierr;
//	PetscInt i, n = 10, col[3], its;
//	PetscMPIInt size;
//	PetscScalar neg_one = -1.0, one = 1.0, value[3];
//	PetscBool nonzeroguess = PETSC_FALSE;

	//ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
	//if (size != 1) SETERRQ(PETSC_COMM_WORLD, 1, "This is a uniprocessor example only!");
	//ierr = PetscOptionsGetInt(PETSC_NULL, "-n", &n, PETSC_NULL); CHKERRQ(ierr);
	//ierr = PetscOptionsGetBool(PETSC_NULL, "-nonzero_guess", &nonzeroguess, PETSC_NULL); CHKERRQ(ierr);

	//ierr = MatCreate(PETSC_COMM_WORLD, &A);



	//This can be not necessary
	//vector<Element> ve; ve.push_back(e);
	//FeGrid grid(ve);




	logfile.close();
	return 0;


}