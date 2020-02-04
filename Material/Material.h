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

#ifndef _MATERIAL_H_
#define _MATERIAL_H_

namespace FluxSol
{


	class Material{

	
	protected:
	unsigned int id;
	
	vector <double> k;	//Constants. TO MODIFY, CAN BE FUNCTIONS
	public:
	Material(){}		//Generic Constant Material Constructor
	
	Material (int klen, ...)
	{
		va_list vl;
		va_start(vl, klen);

		for (int i = 0; i < klen; ++i)
		{
				k.push_back(va_arg(vl, double));
		}
		va_end(vl);
	}
	
	virtual ~Material(){}

	const double & K(const int i)const{cout << "K size: "<< k.size()<<endl; return k[i];}
	};



	class LinearElasticMaterial :
		public Material
	{
		

		protected:

		public:
			LinearElasticMaterial(){}
			LinearElasticMaterial(int klen, ...)
			{
				va_list vl;
				va_start(vl, klen);

				for (int i = 0; i < klen; ++i)
				{
						k.push_back(va_arg(vl, double));
				}
				va_end(vl);
			}
	};


	//Thermal ISOTROPIC Material
	class ThermalMaterial :
	public Material
	{

		protected:

		public:
			ThermalMaterial(const double &ki) {}

	};


}


#endif