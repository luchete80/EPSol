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
#ifndef _VEC3D_H
#define _VEC3D_H

#include <iostream>
#include <vector>
#include "Scalar.h"
#include "Tensor.h"
#include "Products.h"    //Para incluir tipo de Rango
using namespace std;

//La clase se coloca
// CAMBIAR CARACTERISTICAS DE VEC3D
namespace FluxSol{

//Ambos son del tipo tensor!!
class Vec3D :public Tensor{

private:

protected:
    std::vector<double> angle;   //in radians

public:
	double comp[3];

    enum
    {
        rank = 1 // Rank of Vector is 1, esto sirve para las templates
    };
	//Constructores
	//Vec3D(const double &v){comp[0]=comp[1]=comp[2]=v;};
	Vec3D(double x=0., double y=0., double z=0.)   //Argunentos por defecto
	{
		comp[0]=x;
		comp[1]=y;
		comp[2]=z;
		//Set_angles();
	}
	//Funcion tipica usada en los templates de operaciones
	Vec3D Val();

	//Operadores
	//Is exclusive of vectors
	double dot(const Vec3D &right);
	Vec3D cross(const Vec3D &right);

	//TO MODIFY, linking problems
	Vec3D & operator= (const Vec3D &right)
	{
		comp[0]=right.comp[0];
		comp[1]=right.comp[1];
		comp[2]=right.comp[2];
		return *this;
	};
	Vec3D & operator= (const double &);
	Vec3D & operator= (const std::vector<double> &);
	Vec3D & operator*= (const double &);


	Vec3D &operator/= (const double &);
	Vec3D &operator+= (const double &);
	Vec3D &operator+= (const Vec3D &);
	Vec3D &operator-= (const double &);
	Vec3D &operator-= (const Vec3D &);


	const Vec3D operator* (const double &)const;
	const Vec3D operator/ (const double &)const;
	const Vec3D operator- (const double &)const;
	const Vec3D operator+ (const double &)const;

	const Vec3D operator* (const Vec3D &)const;
	const Vec3D operator/ (const Vec3D &)const;
	const Vec3D operator- (const Vec3D &)const;
	const Vec3D operator+ (const Vec3D &)const;

	const Scalar operator&(const Vec3D &right) const;

	//template<>
	//typename innerProduct<Vec3D, Vec3D> ::type
	//operator& (const Vec3D &right);

	bool operator== (const Vec3D &);
	bool operator!= (const Vec3D &);
	double & operator[] (int i){return comp[i];};
	const Scalar Norm()const;							//This must be moved to operations

	const Vec3D normalize(void) ;	//Versor


    //Componentes
    inline const double & x()const{return comp[0];}
    inline const double & y()const{return comp[1];}
    inline const double & z()const{return comp[2];}


	vector<double> Comp();	//Common to all templates

	//Operaciones de Log``
	std::string Imprimir_Coord();


    virtual ~Vec3D(){}; //Destructor virtual
	void Log(ofstream &f){f<<comp[0]<<";"<<comp[1]<<";"<<comp[2]<<";";}

    void Set_angles();
	const double & Angle(const int &i)const{return this->angle[i];}
    const std::vector<double> & Angle()const{return this->angle;}
};

Vec3D  operator^ (const Vec3D &left, const Vec3D &right);

//Funciones aparte
//Puede ser vertice o vec3D
//Vec3D  Prom(const vector <Vec3D> &vecs);
//Para que acepte solo vectores de clase derivadas
//Pero la clase base podria no ser template
template <typename T>
//typename enable_if<is_base_and_derived<BaseClass, T>, unsigned>::type
T  Prom(const vector<T>& vecs)
{
    T vec(0,0,0);
    for (int i=0;i<vecs.size();i++)
        vec+=vecs[i];
    vec/=vecs.size();
    return vec;
}


template<>
class pTraits<Vec3D>
{

    public:
    enum
    {
        //Que dimension es??
        dim = 3,
        rank = 1,
        nComponents = 3
    };

//- Component type
    typedef Vec3D cmptType;
};


//Especializacion de plantilla
template<>
class typeOfRank<1>
{
public:

typedef Vec3D type;
};

const Vec3D operator* (const double &left, const Vec3D &right);
//Inner Prod

//Pending for allocating in Operations.h/cpp
const Vec3D operator* (const Scalar &sc, const Vec3D &vector);
const Vec3D operator* (const Vec3D &vector, const Scalar &sc);

Vec3D operator/ (Vec3D &vector, Scalar &sc);
double fabs(const Vec3D& vec);

}//Fin de namespace FluxSol

//Operaciones que no son de clase pero que incluyen a otras clases



//Funciones inline
#include "Vec3dI.h"

#endif
