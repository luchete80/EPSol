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
namespace FluxSol{
//Puedo colocarlo como miembro
//PASARLAS A INLINE!!

//double fabs(const Vec3D& vec);
// This operators must not be
const Vec3D operator/ (const double &left, const Vec3D &right);
const Vec3D operator+ (const double &left, const Vec3D &right);
const Vec3D operator- (const double &left, const Vec3D &right);


ostream &operator<< (ostream &output,const Vec3D &right);

}//Fin FluxSol
