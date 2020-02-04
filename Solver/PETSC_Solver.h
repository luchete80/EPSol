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

#ifndef _PETSC_SOLVER_H_
#define _PETSC_SOLVER_H_

#include "Solver.h"

#include "petscksp.h"
#include <vector>

#include "Vector.h"
#include "DoFHandler.h"

using namespace std;

namespace FluxSol
{

//TO MODIFY: MUST INCLUDE A EQN SYSTEM
//THIS ONTE CAN BE TEMPLATIZeD WITH MAtrIX; SPARSE MATRIX, MAT
//INHERIT FROM THIS SOLVER, SEQUENTIAL AND PARALLEL SOLVERS
template <typename number, int dim>
class PETSC_Solver:
public Solver<number, dim>
{

	protected:

		// PETSC variables
		Mat A;
		Vec b; // right hand side
		Vec x; // solution, residual vectors
		Mat SysMat;

		KSP ksp; // linear solver context
		PC pc; // preconditioner context

		PetscErrorCode ierr;

		PetscMPIInt size;
		PetscScalar neg_one = -1.0,one = 1.0;
		PetscBool nonzeroguess = PETSC_FALSE;

        //TO MODIFY
		//CAN HAVE A DOFHANDLER
		//DoFHandler<dim> &dofhandler;

		//TO MODIFY
		//REFERENCE TO DOFHANDLER OR GRID!!

	public:

	void PETSC_Init();
	//Constructors
	PETSC_Solver():Solver<number,dim>(0){};
	PETSC_Solver(const int &d);

	void PreAllocateRows(const vector<int> &nnz);
	void PreAllocateRows(const PetscInt &cols);

	void Solve();
	void InsertRow(const int &row, const std::vector<int> &cols, const std::vector <double> &vals);
	void ResetMatrix(const DoFHandler<dim> &dofhandler);
	void ResetB(const DoFHandler<dim> &dofhandler);

	inline void SetMatVal(const PetscInt &row, const PetscInt &col, const PetscScalar &value);
	inline void AddMatVal(const PetscInt &row, const PetscInt &col, const PetscScalar &value);
	void SetbValues(const PetscInt &row, const PetscScalar &value);	//TO MODIFY, TEMPLATIZE TO VALUE
	void AddbValues(const PetscInt &row, const PetscScalar &value);	//TO MODIFY, TEMPLATIZE TO VALUE
	void SetbValues(const PetscScalar &value);
	void ApplyBCOnDoF(const int &dof, const DoFHandler<dim> &dofhandler);
    void ApplyBCOnDoF(const int &dof, const DoFHandler<dim> &dofhandler, const number &n);
	void ApplyDispOnDoF(const int &dof, const number &u, const DoFHandler<dim> &dofhandler);

	const FluxSol::Vector <number> X() const;	//Returns X Solution
	const number X(const int &i) const;	//Returns X Solution
	const FluxSol::Vector <number> B() const;
	const number B(const int &i) const;

	void ViewInfo();

	void Flush(){
		ierr = MatAssemblyBegin(this->A,MAT_FLUSH_ASSEMBLY);
		ierr = MatAssemblyEnd  (this->A,MAT_FLUSH_ASSEMBLY);
	}

	inline void ClearMat();

	// ierr = VecDestroy(&x);CHKERRQ(ierr); ierr = VecDestroy(&u);CHKERRQ(ierr);
	// ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
	// ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
	// /*
	// Always call PetscFinalize() before exiting a program. This routine
	// - finalizes the PETSc libraries as well as MPI
	// - provides summary and diagnostic information if certain runtime
	// options are chosen (e.g., -log_summary).
	// */
	// ierr = PetscFinalize();

};//PETSC Solver

}//FluxSol
#endif
