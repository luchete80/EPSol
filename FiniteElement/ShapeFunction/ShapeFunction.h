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

#ifndef _SHAPEFUNCTION_H_
#define _SHAPEFUNCTION_H_

#include "Scalar.h"
#include "../../Function/Function.h"
#include "Integration\GaussPoint.h"

#include <sstream>      // std::ostringstream
#include <iomanip>      // std::setprecision
#include <cstdarg>

#include <string>       // std::string
#include <iostream>     // std::cout

#include <vector>

using namespace std;

namespace FluxSol
{

	//Explicit template
	//IORIGINALLY
	//IF IT IS QUADRATIC CAN INHERIT MEMBER VARIABLES
	//ACTUALLY CAN DIFERENCIATE BETWEEN THEM
	class ShapeFunction :
		public PropFunction<double>
	{

	protected:


	public:
		ShapeFunction(){}

		//linear shape function
		ShapeFunction(const int &order) :PropFunction(order){}
		ShapeFunction(	const double &r, const double &s, const double &t,
						const double &rs, const double &st, const double &rt,
						const double &kt) :PropFunction(2)
						{
							this->coeff[1][0][0] = r; this->coeff[0][1][0] = s; this->coeff[0][0][1] = t;
							this->coeff[1][1][0] = rs; this->coeff[0][1][1] = st; this->coeff[1][0][1] = rt;
							this->coeff[0][0][0] = kt;
						}

		//Same as previous
		//Added coeff rst
        ShapeFunction(	const double &r, const double &s, const double &t,
                        const double &rs, const double &st, const double &rt,
                        const double &kt,
                        const double &rst) :
                        PropFunction(2)
						{
						    this->coeff[1][0][0] = r; this->coeff[0][1][0] = s; this->coeff[0][0][1] = t;
							this->coeff[1][1][0] = rs; this->coeff[0][1][1] = st; this->coeff[1][0][1] = rt;
							this->coeff[0][0][0] = kt;
							this->coeff[1][1][1] = rst;
						}
		//Can Make an overload operator[][][] Like in Matrix
		//NO CONST
		//TO MODIFY. GP Function Vals may be storaged
		inline virtual const double Val(const double &r, const double &s, const double &t) const;
		//const double & r()
		//Evaluate in gausspoint
		const double Val(const GaussPoint &gp)const { return Val(gp.r(),gp.s(),gp.t()); };
		//DIFFERS FROM LINEAR TO CUBIC FUNCTIONS
		// NON VIRTUAL!!!
		inline const vector<ShapeFunction> diff() const;	//returns a vector of inherited functions

		inline std::string outstr() const;
	};



	//Explicit template
	//THIS MUST INHERIT FROM LINEAR FUNCTION
	class LinearShapeFunction :
		public ShapeFunction
	{


	public:
		LinearShapeFunction(){}

		inline const double Val(const double &r, const double &s, const double &t) const;


	};


	//Explicit template
	class ShapeFunctionGroup
	{
		std::vector<ShapeFunction> shapefn;

	public:
		ShapeFunctionGroup(){}
		ShapeFunctionGroup(const std::vector<ShapeFunction> &shf) :shapefn(shf){}
		void AddShapeFn(const ShapeFunction &f){ shapefn.push_back(f); };

		//Overloading operator to return single function
		ShapeFunction & operator[](const int &i){ return shapefn[i]; }
		const double Size()const { return shapefn.size(); }
		inline const ShapeFunction& ShapeFn(const unsigned int &i)const{ return shapefn[i]; }

		//DIFF IS NOT VIRTUAL
		//const virtual vector<ShapeFunctionGroup> diff()const;	//returns a vector of inherited functions
	};


	inline const double ShapeFunction::Val(const double &r, const double &s, const double &t)const
	{
		double res;
		//T kterm;		//Constant term
		//std::vector <T> linearcoeff;	//Six terms, r,s,t and cross

		res = r*coeff.Val(1,0,0) + s*coeff.Val(0,1,0) + t *coeff.Val(0,0,1) +
			r*s*coeff.Val(1,1,0) + s*t*coeff.Val(0,1,1) + r*t*coeff.Val(1,0,1) +
			coeff.Val(0,0,0)+
			r*s*t*coeff.Val(1,1,1);

		const double ret = res;

		return ret;

	}

	// NONVIRTUAL
	const vector<ShapeFunction> ShapeFunction::diff() const
	{
		ShapeFunction fr(this->order_);
		ShapeFunction fs(this->order_);
		ShapeFunction ft(this->order_);


		ShapeFunction temp = *this;
		//If exponent>=0, derivate
		for (int var1 = 0; var1 < order_; var1++)
		{
			for (int var2 = 0; var2 < order_; var2++)
			{
				for (int r = 1; r < order_; r++)
					fr.coeff[r - 1][var1][var2] = r*temp.coeff[r][var1][var2];
				for (int s = 1; s < order_; s++)
					fs.coeff[var1][s - 1][var2] = s*temp.coeff[var1][s][var2];
				for (int t = 1; t < order_; t++)
					ft.coeff[var1][var2][t - 1] = t*temp.coeff[var1][var2][t];
			}
		}
		vector<ShapeFunction> vf;
		vf.push_back(fr); vf.push_back(fs); vf.push_back(ft);


		return vf;
	}

	std::string ShapeFunction::outstr() const
	{

		string st;


		for (int r = 0; r < order_;r ++)
			for (int s= 0; s < order_; s++)
			for (int t = 0; t < order_; t++)
			{
				std::ostringstream strs;
				strs << "r " << r << ",s " << s << ", t " << t << ":  " << this->coeff.Val(r,s,t)<< "\n";
				st+= strs.str();

			}


		return st;
	}

}


#endif
