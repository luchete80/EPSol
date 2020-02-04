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
#ifndef _INTEGRATIONSCHEME_H_
#define _INTEGRATIONSCHEME_H_


#include "GaussPoint.h"
#include "Matrix.h"
#include "Vector.h"
#include <math.h>
namespace FluxSol{

	template<typename T>
	class IntegrationScheme{

	public:
		inline Matrix<T> Int(const Matrix<T> &m);
		T Int(const T &t);

	};

	class GaussIntegrationScheme
		:public IntegrationScheme<double>
	{
		const int gausspointsnumber;	//CAN NOT BE CHANGED?
		std::vector <GaussPoint> gausspoints;
		int dim;    //TO MODIFY, NEXT TEMPLATIZE

	public:

		//Enum with schemes ids
		//DEFAULT, 2D 4 points //BYNOW
//		GaussIntegrationScheme() :gausspointsnumber(4){
//			vector<GaussPoint> gpv;
//			//Change This
//			gausspoints.push_back(GaussPoint(-0.577350269,-0.577350269, 0.0, 1.0,1.0,0.0));
//			gausspoints.push_back(GaussPoint(-0.577350269, 0.577350269, 0.0, 1.0,1.0,0.0));
//			gausspoints.push_back(GaussPoint( 0.577350269, 0.577350269, 0.0, 1.0,1.0,0.0));
//			gausspoints.push_back(GaussPoint( 0.577350269,-0.577350269, 0.0, 1.0,1.0,0.0));
//		}

        //TO MODIFY
		GaussIntegrationScheme() :gausspointsnumber(0){
		    cout << "ERROR: No dimension set"<<endl;
			vector<GaussPoint> gpv;
			//Change This
			double gpc=0.0;
			double wc =0.888888888888889;
			double gps=0.774596669241483;
			double ws =0.555555555555556;
			//TO MODIFY, DIMENSION ON GAUSS POINTS
			gausspoints.push_back(GaussPoint(-gps,-gps,0.,ws,ws,1.));
			gausspoints.push_back(GaussPoint( gpc,-gps,0.,wc,ws,1.));
			gausspoints.push_back(GaussPoint( gps,-gps,0.,ws,ws,1.));
			gausspoints.push_back(GaussPoint(-gps, gpc,0.,ws,wc,1.));
			gausspoints.push_back(GaussPoint( gpc, gpc,0.,wc,wc,1.));
            gausspoints.push_back(GaussPoint( gps, gpc,0.,ws,wc,1.));
            gausspoints.push_back(GaussPoint(-gps, gps,0.,ws,ws,1.));
			gausspoints.push_back(GaussPoint( gpc, gps,0.,wc,ws,1.));
			gausspoints.push_back(GaussPoint( gps, gps,0.,ws,ws,1.));
		}

		inline GaussIntegrationScheme(const int &order, const int &d):gausspointsnumber(pow((order+1),(double)d))
		{
		    dim=d;

            if (order ==0)
            {
                gausspoints.push_back(GaussPoint(0.,0.,0.,2.,1.,1.));
            }
            else if (order ==1)
            {
                double p=1.0/1.732050807568877;
                gausspoints.push_back(GaussPoint(-p,-p,0.,1.,1.,1.));
                gausspoints.push_back(GaussPoint( p,-p,0.,1.,1.,1.));
                gausspoints.push_back(GaussPoint(-p, p,0.,1.,1.,1.));
                gausspoints.push_back(GaussPoint( p, p,0.,1.,1.,1.));
            }
            else if (order == 2)//3 Gauss Points
            {


                Vector<double> gpgen(3);
                Vector<double> wgen(3);

                gpgen[0]=-0.774596669241483;
                gpgen[1]= 0.;
                gpgen[2]= 0.774596669241483;

                wgen[0]= 0.555555555555556;
                wgen[1]= 0.888888888888889;
                wgen[2]= 0.555555555555556;
                Vec3D vc;
                Vec3D vw;




                //TO MODIFY
                if (dim==1)
                {

                    for (int r=0;r<3;r++)
                    {
                        vc=Vec3D(gpgen[r],0.,0.);
                        vw=Vec3D(wgen[r],0.,0.);
                        gausspoints.push_back(GaussPoint(vc,vw));
                    }


                }
                else if (dim==2)
                {

                    for (int s=0;s<3;s++)
                        for (int r=0;r<3;r++)
                        {
                            vc=Vec3D(gpgen[r],gpgen[s],0.);
                            vw=Vec3D(wgen[r],wgen[s],1.);
                            gausspoints.push_back(GaussPoint(vc,vw));
                        }
                }
                else if (dim==3)
                {

                    for (int t=0;t<3;t++)
                        for (int s=0;s<3;s++)
                            for (int r=0;r<3;r++)
                            {
                                vc=Vec3D(gpgen[r],gpgen[s],gpgen[t]);
                                vw=Vec3D(wgen[r],wgen[s],wgen[t]);
                                gausspoints.push_back(GaussPoint(vc,vw));
                            }

                }


            }//if order 2

		}
		//DIMENSION 1 TO 3
		//inline GaussIntegrationScheme(const int &order, const int &d) {};

		//inline FESparseMatrix Int(const FESparseMatrix m);
		const GaussPoint & operator[](const int &i)const{ return gausspoints[i]; }

		inline const int & NumPoints() const { return this->gausspointsnumber; }

		inline GaussIntegrationScheme & operator=(const GaussIntegrationScheme &is){
			for (int g = 0; g < is.NumPoints(); g++)
				this->gausspoints[g] = is.gausspoints[g];
            this->dim=is.dim;
			//this->gausspointsnumber = is.gausspointsnumber;
			return *this;
		}

	};

}

#endif
