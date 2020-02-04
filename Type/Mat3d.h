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
#ifndef _MAT3D_H
#define _MAT3D_H


#include <vector>
#include <fstream>
#include "../Type/Tensor.h"
#include "../Type/pTraits.h"
#include "../Type/Products.h"


namespace FluxSol
{
	class Mat3D :		//Matrix like Tensor  3 x 3
		public Tensor

	{

		double v_[9];

	public:

		enum
		{
			rank = 2 // Rank of Vector is 1, esto sirve para las templates
		};
		//Like Tensor in OpenFoam
		enum components { XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ };

		Mat3D(){}
		inline Mat3D
			(
			const double txx, const double txy, const double txz,
			const double tyx, const double tyy, const double tyz,
			const double tzx, const double tzy, const double tzz
			);

		//Constructors


		inline const double& xx() const;
		inline const double& xy() const;
		inline const double& xz() const;
		inline const double& yx() const;
		inline const double& yy() const;
		inline const double& yz() const;
		inline const double& zx() const;
		inline const double& zy() const;
		inline const double& zz() const;

		inline double& xx();
		inline double& xy();
		inline double& xz();
		inline double& yx();
		inline double& yy();
		inline double& yz();
		inline double& zx();
		inline double& zy();
		inline double& zz();

		//Member operators
		inline Mat3D operator/(const double &d);

	};



	//Original in OpenFOam
	//template<class Cmpt>
	//class typeOfRank<Cmpt, 2>
	template<>
	class typeOfRank<2>
	{
	public:
		typedef Mat3D type;

	};


	//Inline
	//Constructors
	//- Construct from components
	inline Mat3D::Mat3D
		(
		const double txx, const double txy, const double txz,
		const double tyx, const double tyy, const double tyz,
		const double tzx, const double tzy, const double tzz
		)
	{
		this->v_[XX] = txx; this->v_[XY] = txy; this->v_[XZ] = txz;
		this->v_[YX] = tyx; this->v_[YY] = tyy; this->v_[YZ] = tyz;
		this->v_[ZX] = tzx; this->v_[ZY] = tzy; this->v_[ZZ] = tzz;
	}


	//Member operators

	inline Mat3D Mat3D::operator/(const double &d)
	{


	}



	//- Return the determinant of a tensor
	//TO MODIFY: GENERALIZE
	inline double det(const Mat3D& t)
	{
		return
			(
			t.xx()*t.yy()*t.zz() + t.xy()*t.yz()*t.zx()
			+ t.xz()*t.yx()*t.zy() - t.xx()*t.yz()*t.zy()
			- t.xy()*t.yx()*t.zz() - t.xz()*t.yy()*t.zx()
			);
	}


	//- Return the cofactor tensor of a tensor
	inline Mat3D cof(const Mat3D& t)
	{
		return Mat3D
			(
			t.yy()*t.zz() - t.zy()*t.yz(),
			t.zx()*t.yz() - t.yx()*t.zz(),
			t.yx()*t.zy() - t.yy()*t.zx(),

			t.xz()*t.zy() - t.xy()*t.zz(),
			t.xx()*t.zz() - t.xz()*t.zx(),
			t.xy()*t.zx() - t.xx()*t.zy(),

			t.xy()*t.yz() - t.xz()*t.yy(),
			t.yx()*t.xz() - t.xx()*t.yz(),
			t.xx()*t.yy() - t.yx()*t.xy()
			);
	}


	//- Return the inverse of a tensor given the determinant
	inline Mat3D inv(const Mat3D& t, const double dett)
	{
		return Mat3D
			(
			t.yy()*t.zz() - t.zy()*t.yz(),
			t.xz()*t.zy() - t.xy()*t.zz(),
			t.xy()*t.yz() - t.xz()*t.yy(),

			t.zx()*t.yz() - t.yx()*t.zz(),
			t.xx()*t.zz() - t.xz()*t.zx(),
			t.yx()*t.xz() - t.xx()*t.yz(),

			t.yx()*t.zy() - t.yy()*t.zx(),
			t.xy()*t.zx() - t.xx()*t.zy(),
			t.xx()*t.yy() - t.yx()*t.xy()
			) / dett;
	}


	//- Return the inverse of a tensor
	inline Mat3D inv(const Mat3D& t)
	{
		return inv(t, det(t));
	}

}//FluxSol

#endif
