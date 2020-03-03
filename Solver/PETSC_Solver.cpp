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

#include "PETSC_Solver.h"

namespace FluxSol

{

template <typename number,int dim>
void PETSC_Solver<number,dim>::PETSC_Init()
{
	int n=this->matdim;

	cout << "Initializing Solver..."<<endl;

	//If no values passed
	int argc;
	char **args;

	char help[100];

	argc=0;
	args=NULL;

	PetscInitialize(&argc,&args,(char *)0,help);


	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
	if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");
	ierr = PetscOptionsGetInt(PETSC_NULL,"-n",&n,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(PETSC_NULL,"-nonzero_guess",&nonzeroguess,PETSC_NULL);CHKERRQ(ierr);

	/*
	Create vectors. Note that we form 1 vector from scratch and
	then duplicate as needed.
	*/
	cout <<"Creating matrix and vectors"<<endl;
	ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
	ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
	ierr = VecSetFromOptions(x);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
	//ierr = VecDuplicate(x,&u);CHKERRQ(ierr);

	/*
	Create matrix. When using MatCreate(), the matrix format can
	be specified at runtime.
	Performance tuning note: For problems of substantial size,
	preallocation of matrix memory is crucial for attaining good
	performance. See the matrix chapter of the users manual for details.
	*/

	ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);	//Instead of create Sij
	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
	ierr = MatSetFromOptions(A);CHKERRQ(ierr);
	//TEMPORARY
	//THIS IS NOT RECOMMENDED DIRECTLY
	//ierr = MatCreateSeqBAIJ(PETSC_COMM_WORLD,PetscInt bs,PetscInt m,PetscInt n,PetscInt nz,const PetscInt nnz[],Mat *A)
	//SEQSBAIJ

	//Symmetric, THIS IS TEMP, TO MODIFY
	//ierr = MatSetOption(A,MAT_SYMMETRIC,PETSC_TRUE);

	//Must call MatXXXSetPreallocation() or MatSetUp() on argument 1 "mat" before MatSetValues()!

	ierr = MatSetUp(A);
	//
	//ierr= MatSeqAIJSetPreallocation(Mat B,PetscInt nz,const PetscInt nnz[])
	//where
	//B	- The matrix-free
	//nz	- number of nonzeros per row (same for all rows)
	//nnz	- array containing the number of nonzeros in the various rows (possibly different for each row) or NULL



	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Create the linear solver and set various options
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/*
	Create linear solver context
	*/
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	/*
	Set operators. Here the matrix that defines the linear system
	also serves as the preconditioning matrix.
	*/

	ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
	/*
	Set linear solver defaults for this problem (optional).
	- By extracting the KSP and PC contexts from the KSP context,
	we can then directly call any KSP and PC routines to set
	various options.
	- The following four statements are optional; all of these
	parameters could alternatively be specified at runtime via
	KSPSetFromOptions();
	*/
	ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
	ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
	//KSPSetTolerances(KSP ksp,PetscReal rtol,PetscReal abstol,PetscReal dtol,PetscInt maxits)
	ierr = KSPSetTolerances(ksp,1.e-50,1.e-50,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	/*
	Set runtime options, e.g.,
	-ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
	These options will override those specified above as long as
	KSPSetFromOptions() is called _after_ any other customization
	routines.
	*/
	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
	if (nonzeroguess)
	{
		PetscScalar p = .5;
		ierr = VecSet(x,p);CHKERRQ(ierr);
		ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);
	}

	ierr=VecSet(this->b,0.);
	ierr=VecSet(this->x,0.);
	cout << "Solver Initialized." <<endl;
}

template <typename number,int dim>
PETSC_Solver<number,dim>::
PETSC_Solver(const int &d):
Solver<number,dim>(d)
{
	PETSC_Init();
}


template <typename number,int dim>
void PETSC_Solver<number,dim>::PreAllocateRows(const vector <int> &nnz)
{
	MatSeqAIJSetPreallocation(this->A,PETSC_NULL,&nnz[0]);
}

template <typename number,int dim>
void PETSC_Solver<number,dim>::PreAllocateRows(const PetscInt &cols)
{
	MatSeqAIJSetPreallocation(this->A,cols,PETSC_NULL);
}

template <typename number,int dim>
void PETSC_Solver<number,dim>::Solve()
{
	ierr = MatAssemblyBegin(this->A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(this->A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	//ierr = MatSetOption(this->A,MAT_SYMMETRIC,PETSC_TRUE);


	ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
	/*
	View solver info; we could instead use the option -ksp_view to
	print this info to the screen at the coknclusion of KSPSolve().
	*/
	ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Check solution and clean up
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/*
	Check the error
	*/
	PetscReal norm;
	PetscInt its;

	ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
	ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %A, Iterations %D\n",
	norm,its);CHKERRQ(ierr);
	/*
	Free work space. All PETSc objects should be destroyed when they
	are no longer needed.
	*/

}

template <typename number,int dim>
void PETSC_Solver<number,dim>::InsertRow(const int &row, const std::vector<int> &cols, const std::vector <double> &vals)
{
	//ierr = MatSetValues(this->A,1,row,cols.size(),cols,&vals[0],INSERT_VALUES);CHKERRQ(ierr);
}

template <typename number,int dim>
void PETSC_Solver<number,dim>::SetMatVal(const PetscInt &row, const PetscInt &col, const PetscScalar &value)
{
	ierr=MatSetValues(this->A,1,&row,1,&col,&value,INSERT_VALUES);
}

template <typename number,int dim>
void PETSC_Solver<number,dim>::AddMatVal(const PetscInt &row, const PetscInt &col, const PetscScalar &value)
{
	ierr=MatSetValues(this->A,1,&row,1,&col,&value,ADD_VALUES);
}


template <typename number,int dim>
void PETSC_Solver<number,dim>::ViewInfo()
{

	ierr = MatAssemblyBegin(this->A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(this->A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	cout << "Matrix A: "<<endl;
	MatView(this->A,PETSC_VIEWER_STDOUT_SELF);
	cout << "RHS: "<<endl;
	VecView(b,PETSC_VIEWER_STDOUT_SELF);
	cout << "Solution: "<<endl;
	VecView(x,PETSC_VIEWER_STDOUT_SELF);
}

template <typename number,int dim>
void PETSC_Solver<number,dim>::SetbValues(const PetscInt &row, const PetscScalar &value)
{
	ierr=VecSetValues(this->b,1,&row,&value,INSERT_VALUES);
}

template <typename number,int dim>
void PETSC_Solver<number,dim>::AddbValues(const PetscInt &row, const PetscScalar &value)
{
	ierr=VecSetValues(this->b,1,&row,&value,ADD_VALUES);
}

template <typename number,int dim>
void PETSC_Solver<number,dim>::SetbValues(const PetscScalar &value)
{
	ierr=VecSet(this->b,value);
}

//This function must not call get_dof_indices, because is called directly dof
//Applies a fixity
template <typename number,int dim>
void PETSC_Solver<number,dim>::ApplyBCOnDoF(const int &dof, const DoFHandler<dim> &dofhandler)
{

	for (int i = 0; i < dofhandler.Adj_DoF_Number()[dof]; i++)
	{
		unsigned int adof = dofhandler.AdjDoF(dof,i);
		//cout << "Applying BC on DoF [dof][i]: "<< dof << " "<<i <<": "<< adof<<endl;
		//Have to set up only the nonzero matrix values
		this->SetMatVal(dof,adof,0.);
		this->SetMatVal(adof,dof,0.);
	}
	this->SetMatVal(dof,dof,1.);
	this->SetbValues(dof, 0.);


}

template <typename number,int dim>
void PETSC_Solver<number,dim>::ApplyBCOnDoF(const int &dof, const DoFHandler<dim> &dofhandler, const number &n)
{
    number bval,matval;
	double y[1];
	int idxm[1];    //row
	int idxn[1];    //col

    //Look through all adjacent node numbers
	for (int i = 0; i < dofhandler.Adj_DoF_Number()[dof]; i++)
	{
	    ierr = MatAssemblyBegin(this->A,MAT_FINAL_ASSEMBLY);
        ierr = MatAssemblyEnd  (this->A,MAT_FINAL_ASSEMBLY);

		unsigned int adof = dofhandler.AdjDoF(dof,i);
		//cout << "Applying BC on DoF [dof][i]: "<< dof << " "<<i <<": "<< adof<<endl;
		//Have to set up only the nonzero matrix values
		//bval=this->SysMat
		//ierr=MatGetValues(Mat mat,PetscInt m,const PetscInt idxm[],PetscInt n,const PetscInt idxn[],PetscScalar v[])
        //cout << "Getting values" <<endl;
        idxm[0]=adof;idxn[0]=dof;
        //cout << "Adjdof " <<adof << " Dof " <<dof<<endl;
        ierr=MatGetValues(this->A,1,idxm,1,idxn,y);
        //ierr=MatGetValues(this->A,1,&adof,1,&dof,y);
        matval=y[0];
        //cout << "value "<<matval<<endl;
        bval=-matval*n;
        //cout << "bval "<<bval<<endl;

        ierr = MatAssemblyBegin(this->A,MAT_FLUSH_ASSEMBLY);
        ierr = MatAssemblyEnd  (this->A,MAT_FLUSH_ASSEMBLY);

		this->SetMatVal(dof,adof,0.);
		this->SetMatVal(adof,dof,0.);
        //VecSetValues(Vec x,PetscInt ni,const PetscInt ix[],const PetscScalar y[],InsertMode iora)
        //cout << "Setting b values adj pos "<< adof<<endl;
        ierr=VecSetValues(this->b,1,&adof,&bval,ADD_VALUES);
	}
	this->SetMatVal(dof,dof,1.);
	ierr=VecSetValues(this->b,1,&dof,&n,INSERT_VALUES);



}

template <typename number,int dim>
void PETSC_Solver<number,dim>::ApplyDispOnDoF(const int &dof, const number &u, const DoFHandler<dim> &dofhandler)
{
	ApplyBCOnDoF(dof, dofhandler);
	ierr=VecSetValues(this->b,1,&dof,&u,INSERT_VALUES);

}

template <typename number, int dim>
const FluxSol::Vector <number>
PETSC_Solver<number,dim>::X() const
{
	Vector <number> v(this->matdim);

	int ix[1];
	double y[1];

	for (int i=0;i<this->matdim;i++) //TO EVALUATE SPEED
	{
		ix[0]=i;
		VecGetValues(x,1,ix,y);
		number val;
		val=y[0];
		v[i]=val;
	}

	//PetscErrorCode  VecGetValues(Vec x,PetscInt ni,const PetscInt ix[],PetscScalar y[])

	return v;
}

template <typename number,int dim>
void PETSC_Solver<number,dim>::ClearMat()
{
	//TO MODIFY
	int row,col;
	//for (int i=0;i<)
	//ierr=MatSetValue(this->A,1,&row,1,&col,&value,INSERT_VALUES);

}

template <typename number,int dim>
const number
PETSC_Solver<number,dim>::B(const int &i)const
{
	double y[1];
	int ix[1];
	ix[0]=i;

	VecGetValues(this->b,1,ix,y);
	//cout << "value get" << y[0] <<endl;
	number n=y[0];
	return n;
}

template <typename number,int dim>
const number
PETSC_Solver<number,dim>::X(const int &i)const
{
	double y[1];
	int ix[1];
	ix[0]=i;

	VecGetValues(this->x,1,ix,y);
	//cout << "value get" << y[0] <<endl;
	number n=y[0];
	return n;
}

template <typename number, int dim>
const FluxSol::Vector <number>
PETSC_Solver<number,dim>::B() const
{
	Vector <number> v(this->matdim);

	int ix[1];
	double y[1];

	for (int i=0;i<this->matdim;i++) //TO EVALUATE SPEED
	{
		ix[0]=i;
		VecGetValues(this->b,1,ix,y);
		number val;
		val=y[0];
		v[i]=val;
	}

	//PetscErrorCode  VecGetValues(Vec x,PetscInt ni,const PetscInt ix[],PetscScalar y[])

	return v;
}

template <typename number,int dim>
void PETSC_Solver<number,dim>::ResetMatrix(const DoFHandler<dim> &dofhandler)
{
    for (int dof=0;dof<dofhandler.NumDoF();dof++)
    //int ix[1];
    //TO MODIFY: MAKE IT FAST, CALL PETSC Function MatSetValues once per Matrix Solver row
	for (int i = 0; i < dofhandler.Adj_DoF_Number()[dof]; i++)
	{
		unsigned int adof = dofhandler.AdjDoF(dof,i);
		//cout << "Applying BC on DoF [dof][i]: "<< dof << " "<<i <<": "<< adof<<endl;
		//Have to set up only the nonzero matrix values
		this->SetMatVal(dof,adof,0.);
		this->SetMatVal(adof,dof,0.);
	}
}

template <typename number,int dim>
void PETSC_Solver<number,dim>::ResetB(const DoFHandler<dim> &dofhandler)
{
    for (int dof=0;dof<dofhandler.NumDoF();dof++)
    {
		this->SetbValues(dof,0.);
	}
}

}//FluxSol

#include "PETSC_Solver.inst"
