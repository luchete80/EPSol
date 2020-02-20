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
// EXAMPLE 11 //
///////////////
// Surface Loads Examples
/////////////////////////

#include "EPSol.h"

//Sparse libraries
#include "laspack.h"

#include <iostream>
#include <vector>

using namespace std;
using namespace FluxSol;

int main()
{

	ofstream logfile;
	logfile.open("Logfile.txt");

	//Element can be construted from vertex too
	//std::_Vertex

	std::vector<Node> v;

	//TO MODIFY: "ADDNODEFUNCTION IN GRID"
	v.push_back(Node(0, 1.0, 1.0, 0.0));
	v.push_back(Node(1, 0.0, 1.0, 0.0));
	v.push_back(Node(2, 0.0, 0.0, 0.0));
	v.push_back(Node(3, 1.0, 0.0, 0.0));

	const FEIntegrationScheme intsch;


	QuadLinearElement e(v);
	e.Set_Nodes(4, 0, 1, 2, 3);
	logfile << "Element Gauss Order: " << e.GaussOrder() << "\n\n";

	GaussMatrices H = e.H();

	std::vector< Element<2> > ve;
	ve.push_back(e);
	FeGrid <2> g(v, ve);
	FEValues<2> fev(e, g);

	return 0;


}
