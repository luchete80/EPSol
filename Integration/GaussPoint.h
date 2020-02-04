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

#ifndef _GAUSSPOINT_H_
#define _GAUSSPOINT_H_

#include "../Type/Vec3D.h"
//#include "../Type/Point.h"

namespace FluxSol{

    //TO MODIFY: POINT MUST INHERIT FROM TENSOR
	class GaussPoint
	//:public Point<dim>
	{

        double weight;
		const Vec3D wv;
		const Vec3D coords;
		int order;

        //TO MODIFY
        //TEMPLATIZE!!
		public:
			GaussPoint(){}
			GaussPoint(const Vec3D v, const double &w) :coords (v)       {  weight = w; }
			GaussPoint(const Vec3D &v, const Vec3D &w) :coords (v)
			{
                weight= w.comp[0]*
                w.comp[1]*
                w.comp[2];
            }
			GaussPoint(	const double &x,
						const double &y,
						const double &z,
						const double &wx,
						const double &wy,
						const double &wz):coords(Vec3D(x, y, z)),wv(Vec3D(wx,wy,wz))
						{
                            weight= wv.comp[0]*
                                    wv.comp[1]*
                                    wv.comp[2];

						}

			const double & r()const{ return coords.comp[0]; };
			const double & s()const{ return coords.comp[1]; };
			const double & t()const{ return coords.comp[2]; };

            const Vec3D & Coords()const{return this->coords;}

			inline GaussPoint& operator=(const GaussPoint &g)
			{
				//const Vec3D (g.coords);
				//const Vec3D (g.weight);
				this->order = g.order;
                this->weight = g.order;
				return *this;

			}

			const double & w()const { return this->weight; };
	};


};
#endif
