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
//EXACT LIKE EXAMPLE 1 FROM PETSC, BUT ADDED MatS
/*
Include "petscksp.h" so that we can use KSP solvers. Note that this file
automatically includes:
petscsys.h - base PETSc routines petscvec.h - vectors
petscmat.h - matrices
petscis.h - index sets petscksp.h - Krylov subspace methods
petscviewer.h - viewers petscpc.h - preconditioners
Note: The corresponding parallel example is ex23.c
*/
#include <petscksp.h>
//#undef __FUNCT__
//#define __FUNCT__ "main"
int main(int argc,char **args)
{
	Vec x, b, u; /* approx solution, RHS, exact solution */
	Mat A; /* linear system matrix */
	KSP ksp; /* linear solver context */
	PC pc; /* preconditioner context */
	PetscReal norm; /* norm of solution error */
	PetscErrorCode ierr;

	char help[100];
	
	PetscInt i,n = 3,col[3],its;
	PetscMPIInt size;
	PetscScalar neg_one = -1.0,one = 1.0,value[3];
	PetscBool nonzeroguess = PETSC_FALSE;
	PetscInitialize(&argc,&args,(char *)0,help);
	
	
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
	if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");
	ierr = PetscOptionsGetInt(PETSC_NULL,"-n",&n,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(PETSC_NULL,"-nonzero_guess",&nonzeroguess,PETSC_NULL);CHKERRQ(ierr);
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	24
	Compute the matrix and right-hand-side vector that define
	the linear system, Ax = b.
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/*
	Create vectors. Note that we form 1 vector from scratch and
	then duplicate as needed.
	*/
	ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
	ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
	ierr = VecSetFromOptions(x);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&u);CHKERRQ(ierr);
	/*
	Create matrix. When using MatCreate(), the matrix format can
	be specified at runtime.
	Performance tuning note: For problems of substantial size,
	preallocation of matrix memory is crucial for attaining good
	performance. See the matrix chapter of the users manual for details.
	*/
	ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
	ierr = MatSetFromOptions(A);CHKERRQ(ierr);
	
	//MODIFIED BY LUCIANO
	//Must call MatXXXSetPreallocation() or MatSetUp() on argument 1 "mat" before MatSetValues()!
	ierr = MatSetUp(A);
	/*
	Assemble matrix
	*/    
	//"|1  0  2| |x1|   |7|"<<endl;
    //"|0  1  1| |x2| = |5|"<<endl;
    //"|0  0  1| |x3|   |3|"<<endl;
	//With solution:
    //1;1/2;1

	
	//To put several rows at the same time
	//value[0] = -1.0; value[1] = 2.0; value[2] = -1.0;
	//for (i=1; i<n-1; i++) 
	//{
	//	col[0] = i-1; col[1] = i; col[2] = i+1;
	//	ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
	//}
	//i = n - 1; col[0] = n - 2; col[1] = n - 1;
	//ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
	//i = 0; col[0] = 0; col[1] = 1; value[0] = 2.0; value[1] = -1.0;
	//ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
	
	//MatSetValues(Mat mat,PetscInt m,const PetscInt idxm[],PetscInt n,const PetscInt idxn[],const PetscScalar v[],InsertMode addv)
	//v 	- a logically two-dimensional array of values
	//m, idxm 	- the number of rows and their global indices
	//n, idxn 	- the number of columns and their global indices 
	//n IS THE NUMBER OF ROWS TO BE WRITTEN
	
	i = 0; col[0] = 0; col[1] = 1; col[2]=2;
	value[0]=1.;value[1]=0.;value[2]=2.;
	ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
	i = 1; col[0] = 1; col[1] =2;
	value[0]=1.;value[1]=1.;
	ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
	i = 2; col[0] = 2;
	value[0] = 1.0;
	ierr = MatSetValues(A,1,&i,1,col,value,INSERT_VALUES);CHKERRQ(ierr);
	
	
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	/*
	Set exact solution; then compute right-hand-side vector.
	*/
	//ierr = VecSet(u,one);CHKERRQ(ierr);
	//ierr = MatMult(A,u,b);CHKERRQ(ierr);
	//VecSetValues(Vec x,PetscInt ni,const PetscInt ix[],const PetscScalar y[],InsertMode iora)
	//x 	- vector to insert in
	//ni 	- number of elements to add
	//ix 	- indices where to add
	//y 	- array of values 
	col[0]=0;col[1]=1;col[2]=2;
	value[0]=7.;value[1]=5.;value[2]=3.;
	VecSetValues(b,3,col,value,INSERT_VALUES);
	
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
	ierr = KSPSetTolerances(ksp,1.e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
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
	
	//MODIFIED BY LUCIANO
	MatView(A,PETSC_VIEWER_STDOUT_SELF);
	VecView(b,PETSC_VIEWER_STDOUT_SELF);
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Solve the linear system
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/*
	Solve linear system
	*/
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
	ierr = VecAXPY(x,neg_one,u);CHKERRQ(ierr);
	ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
	ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %A, Iterations %D\n",
	norm,its);CHKERRQ(ierr);
	/*
	Free work space. All PETSc objects should be destroyed when they
	are no longer needed.
	*/
	//MODIFIED BY LUCIANO
	VecView(x,PETSC_VIEWER_STDOUT_SELF);
	VecView(b,PETSC_VIEWER_STDOUT_SELF);
	
	ierr = VecDestroy(&x);CHKERRQ(ierr); ierr = VecDestroy(&u);CHKERRQ(ierr);
	ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
	/*
	Always call PetscFinalize() before exiting a program. This routine
	- finalizes the PETSc libraries as well as MPI
	- provides summary and diagnostic information if certain runtime
	options are chosen (e.g., -log_summary).
	*/
	ierr = PetscFinalize();
	return 0;
}