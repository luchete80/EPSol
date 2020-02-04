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
#include "Scalar.h"

#include <cmath>

namespace FluxSol


{

	//Unary
	const Scalar Scalar::operator- () const
	{
		Scalar temp;
		temp.val = -val;
		return temp;
	}

// linking problem
//	const Scalar & Scalar::operator= (const Scalar &right) {
//		val = right.val;
//		return *this;
//	}

	Scalar & Scalar::operator= (const double & right) {

		val = right;
		return *this;
	}

	Scalar & Scalar::operator= (vector<double> &v) //Common to templates for type conversion
	{
		val = v[0];
		return *this;
	}


	Scalar & Scalar::operator-= (double const& right) {

		val -= right;
		return *this;
	}

	Scalar & Scalar::operator+= (double const& right) {

		val += right;
		return *this;
	}

	Scalar & Scalar::operator-= (const Scalar& right)
	{

		val -= right.val;
		return *this;
	}


	Scalar & Scalar::operator+= (const Scalar& right) {

		val += right.val;
		return *this;
	}

	//No side effects
	const Scalar Scalar::operator+ (const Scalar& right) const
	{
		Scalar temp;
		temp = val + right.val;
		return temp;
	}

	const Scalar Scalar::operator- (const Scalar& right) const
	{
		Scalar temp;
		temp = val - right.val;
		return temp;
	}


	const Scalar Scalar::operator+ (const double & right)const
	{
		Scalar temp;
		temp = val + right;
		return temp;
	}


	const Scalar Scalar::operator/ (const Scalar &right) const
	{
		Scalar temp;
		temp = val / right.val;
		return *this;
	}

	const Scalar Scalar::operator/ (const double &right) const
	{
		Scalar temp;
		temp.Val() = val / right;
		return temp;
	}


	const Scalar Scalar::operator* (const Scalar &right) const
	{
		Scalar temp;
		temp = val*right.val;
		return temp;
	}

	const Scalar Scalar::operator* (const double &right) const
	{
		Scalar temp;
		temp = val*right;
		return temp;
	}


	//inline Scalar &Scalar::operator= (const std::vector<double> &right) {
	////	comp[0]=right[0];
	////	comp[1]=right[1];
	////	comp[3]=right[2];
	//	return *this;
	//}


	vector <double> Scalar::Comp()	//Common all templates
	{
		vector <double> r;

		r.push_back(this->val);
		return r;
	}

	const Scalar Scalar::Norm() const
	{

		Scalar temp(fabs(this->val));
		return temp;
	}

	bool Scalar::operator<(const Scalar &right) const
	{
		bool r = false;
		if (val < right.val) r = true;

		return r;
	}

}// Fin de namespace FluxSol
