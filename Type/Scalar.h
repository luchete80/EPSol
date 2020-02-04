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
#ifndef _SCALAR_H
#define _SCALAR_H

#include <vector>
#include <fstream>
#include "../Type/Tensor.h"
#include "../Type/pTraits.h"
#include "../Type/Products.h"

using namespace std;

namespace FluxSol{

//Rank 0 tensor
class Scalar:public Tensor{

protected:
    double val;
public:

    enum
    {
		dim=3,
        rank = 0, // Rank of Vector is 1
		nComponents = 1 
    };

    //Constructores
    //Scalar(){};
    Scalar(const double &value=0.){val=value;}
    double & Val(){return val;} 
	vector <double> Comp ();	//Common all templates

    ~Scalar(){};

    //Operadores - Operators
	//El operador que devuelve es igual al 
	//Estas funciones deben ser inline - These must be inline
	// Operator = can only be a member function

	//Unary
	const Scalar  operator- ()const;

	//Binary
	//Modify & return lvalue
	const Scalar & operator= (const Scalar &right)
	{
		this->val = right.val;
		return *this;
	}
	Scalar & operator= (vector<double> &v) ; //Common to templates to conversion
	Scalar & operator= (const double &right) ;
		
	Scalar & operator+= (const double &right);
	Scalar & operator+= (const Scalar &right);
	Scalar & operator-= (const double &right);
	Scalar & operator-= (const Scalar &right);
	bool operator<(const Scalar &right)const;


	//No effects in object
	const Scalar  operator+ (const double &right) const;
	const Scalar  operator+ (const Scalar& right) const;
	const Scalar  operator- (const Scalar& right) const;
	const Scalar  operator- (const double& right) const;
	const Scalar  operator/ (const double &right) const;
	const Scalar  operator/ (const Scalar &right) const;
	const Scalar  operator* (const Scalar &right) const;
	const Scalar  operator* (const double &right) const;

	const Scalar Norm()const;

	//inline Scalar &operator= (const std::vector <double> &);
	//Non member operators
	friend const Scalar operator+ (const double &left, Scalar &right)
	{
		Scalar r;
		r.Val()=left+right.Val();
		return r;
	}

	friend const Scalar operator* (const double &left, Scalar &right)
	{
		Scalar r;
		r.Val()=left*right.Val();
		return r;
	}

 	friend const Scalar operator/ (const double &left, Scalar &right)
	{
		Scalar r;
		r.Val()=left/right.Val();
		return r;
	}
	friend const Scalar operator- (const double &left, Scalar &right)	
	{
		Scalar r;
		r.Val()=left-right.Val();
		return r;
	}

	void Log(ofstream &f){f<<val;}

};

//Esto es para las dimensiones
//Especializacion de pTraits
//Estilo OpenFoam
template<>
class pTraits<Scalar>
{
    Scalar p_;

    public:
    enum
    {
        dim = 3,
        rank = 0,
        nComponents = 1
    };

//- Component type
    typedef Scalar cmptType;
};


//Esto es lo que devuelve segun el tipo
//Es la explicitacion (especificacion) de template
template<>
class typeOfRank<0>
{
public:

typedef Scalar type;	//Aca es lque devuelve el tipo Scalar si el ragon es 0
};


}//Fin de namespace FluxSol
#endif
