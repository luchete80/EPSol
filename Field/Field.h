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
#ifndef _FIELD_H
#define _FIELD_H

// Parado en directorio Field //
#include "../Type/Scalar.h"
#include "../Type/Vec3d.h"


//#include "../Input/Input.h"

#include <list>

namespace FluxSol{
//Los campos no dependen de la discretizacion
//

//Campo generico
//De este derivan todos los campos
//Este campo no tiene asociada una malla
//Siempre en los campos internos
template<typename T>
class _Field{


    protected:
    std::vector <T> value;
	int order;
	int numberofvals;
	//This must not be able for patchfields since is confusing
	std::vector <int> idgrid_domain;	//Ya sea vertex, surface o volume
										//Son los id de los vertices que forman parte del dominio


    public:
    _Field(){};
    const int Order(){return order;};
	explicit _Field(const int &numval, const double &value =0.);
	const T & Val(const int &i)const{return value[i];};
	void Val (const int &i,const T &v){value[i]=v;}
	const int & Numberofvals()const{return numberofvals;}
	void  Insert_Idgrid_domain(const int &i){idgrid_domain.push_back(i);}
	int & Idgrid_domain(const int &i){return idgrid_domain[i];}
	void Log();
	

	//Operations

	void operator=(const T &val);
	//It is not const since it can be changed
	//Declare both operators to modify or not field
	T& operator[](const int &x){return value[x];}
	const T& operator[](const int &x)const {return value[x];}
	//void operator=(const double &val);

	_Field<T> operator=(const double &val);

	const _Field<T> & operator&(const _Field<T> &right) const;

};//Field

}//FluxSol



#endif
