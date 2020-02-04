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

#ifndef _EPSOL_H
#define _EPSOL_H


///////////////////

///////////////////
// COMMON FILES ///
///////////////////

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <limits>

//////////////////
// library file //
//////////////////
//#include "./Analysis/Analysis.h"

#include "./Boundary/Boundary.h"

//#include "./Debug/Log.h"
//#include "./Debug/Trace.h"

// FIELD
//NOW FIELD INCLUDE FILES ARE SUBDIVIDED
#include "./Field/Field.h"
#include "./Field/GeometricField.h"
#include "./Field/FieldOperations.h"
#include "./Field/Field_Def.h"

//#include "./Field/FvField.h"

//////////////////
//Finite Volume //
//////////////////



///////////////////
//Finite Element //
///////////////////

enum FEAnalysisType{};

#include "./FiniteElement/Mesh/Element.h"
#include "./FiniteElement/Mesh/FEGrid.h"
#include "./FiniteElement/Mesh/CSRGrid.h"

#include "./FiniteElement/ShapeFunction/ShapeFunction.h"

#include "./FiniteElement/FEMAtrix/FEMatrix.h"

#include "./FiniteElement/FEEqnSystem/FEEqnSystem.h"

#include "./FiniteElement/Integration/FEIntegrationScheme.h"
#include "./FiniteElement/Integration/Integration.h"


#include "./FiniteElement/FEValues/FEValues.h"
#include "./FiniteElement/DoFHandler/DoFHandler.h"

///////////////////////////////


//FUNCTIONS
#include "./Function/Function.h"

// Input //
#include "./Input/Input.h"
#include "./Input/SingleInputFile.h"
#include "./Input/Read_Field_Def.h"

// Mesh //
#include "./Mesh/Cell.h"
#include "./Mesh/Face.h"
#include "./Mesh/Grid.h"
#include "./Mesh/Node.h"
#include "./Mesh/Structured.h"
#include "./Mesh/Vertex.h"

//#include "./Model/Model.h"



// Nastran
//#include "./Nastran/Cadenas.h"
//#include "./Nastran/Nastran.h"
//#include "./Nastran/SistCoord.h"
//#include "./Nastran/Varios.h"



// SistEcuac //
#include "./SistEcuac/SistEcuac.h"
#include "./SistEcuac/SistEcuacDef.h"



//#include "./Tmp/Tmp.h"
//#include "./Tmp/TmpI.h"

#include "Matrix\Matrix.h" //Matrix Of Any size
#include "Matrix\CoeffMatrix.h" //Matrix Of Any size
#include "Matrix\GaussMatrices.h" //Matrix Of Any size

////////////
// Solver //
////////////
#include "./Solver/Solver.h"
#include "./Solver/PETSC_Solver.h"

//////////
// Type //
//////////
#include "./Type/Vector.h" //Vector Of Any size

#include "./Type/Operations.h"
#include "./Type/Products.h"
#include "./Type/pTraits.h"
#include "./Type/Scalar.h"
#include "./Type/Tensor.h"
#include "./Type/Mat3D.h"

//Phisical Tensor
//#include "PTensor.h"

//#include "./Type/TensorI.h"
//Vector
#include "./Type/Vec3d.h"
#include "./Type/Vec3dI.h"

#include "./Type/Mat3d.h"


#include "./Variables/Variable.h"


#endif
