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

#ifndef _VECTOR_H_
#define _VECTOR_H_

#include <vector>
#include <cstdarg>
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::ostringstream
#include <iomanip>      // std::setprecision

using namespace std;

//La clase se coloca
// CAMBIAR CARACTERISTICAS DE Vector
namespace FluxSol{

	//Ambos son del tipo tensor!!
	template <typename T>
	class Vector{

	private:
		vector<T> comp;
		int size;

	public:

		enum
		{
			rank = 1 // Rank of Vector is 1, esto sirve para las templates
		};
		//Constructores
		//Vector(const double &v){comp[0]=comp[1]=comp[2]=v;};
		//EXPLICIT
		inline Vector(int len, ...);   //Argunentos por defecto

		//Funcion tipica usada en los templates de operaciones
		Vector Val();

		const int & Size()const{return this->size;}
		void Clear();

		//Operadores
		//Is exclusive of vectors
		double dot(const Vector &right);
		Vector cross(const Vector &right);
		Vector & operator= (const std::vector<double> &);
		Vector & operator*= (const double &);


		Vector &operator/= (const double &);
		Vector &operator+= (const T &);
		Vector &operator+= (const Vector &);
		Vector &operator-= (const double &);
		Vector &operator-= (const Vector &);


		Vector operator* (const double &);
		Vector operator+ (const double &);

		inline const Vector operator/ (const double &)const;
		inline const Vector operator- (const double &)const;


		inline Vector operator* (const Vector &);
		inline Vector operator/ (const Vector &);
		inline Vector operator- (const Vector &);
		inline Vector operator+ (const Vector &);

		inline const string outstr()const;

		inline Vector<T> & operator=(const Vector<T> &v)
		{
			this->comp.clear();

			for (int i = 0; i < v.Size(); i++)
			{
				comp.push_back(v[i]);
			}
			return *this;
		}


        inline Vector<T> & operator=(const double &d)
		{
			//this->comp.clear();

			for (int i = 0; i < this->Size(); i++)
				comp[i]=d;

			return *this;
		}



		//const Scalar operator&(const Vector &right) const;

		//template<>
		//typename innerProduct<Vector, Vector> ::type
		//operator& (const Vector &right);

		bool operator== (const Vector &);
		bool operator!= (const Vector &);
		inline T & operator[] (int i){return this->comp[i];}
		//const Scalar Norm()const;							//This must be moved to operations

		const Vector normalize(void);	//Versor



		vector<double> Comp();	//Common to all templates

		//Operaciones de Log``
		std::string Imprimir_Coord();


		virtual ~Vector(){}; //Destructor virtual
		void Log(ofstream &f){ f << comp[0] << ";" << comp[1] << ";" << comp[2] << ";"; }

	};




	//Constructor
	template <typename T>
	inline Vector<T>::Vector(int len, ...)
	{
		//std::vector<double> v;
		va_list vl;
		va_start(vl, len);
		comp.push_back(va_arg(vl, T));
		for (int i = 1; i < len; ++i)
			comp.push_back(va_arg(vl, T));
		for (int i = 0; i < len; ++i)
			this->comp[i]=0.;

		this->size=comp.size();
		va_end(vl);
	}


	//Inner product
	// template <typename T>
	// const Vector
	// Vector<T>::operator* (const Vector &)const
	// {
		// Vector<T> ret;

	// }

	template <typename T>
    void Vector<T>::Clear()
    {
        T t=0.;
		for (int c = 0; c < this->comp.size(); c++)
			this->comp[c] = t;
    }

	template <typename T>
	Vector <T> &
	Vector<T>::operator+= (const Vector &v)
	{
		for (int c = 0; c < this->comp.size(); c++)
			this->comp[c] += v[c];

		return *this;
	}


	template <typename T>
	Vector <T>
	Vector<T>::operator* (const double &d)
	{
        Vector <T> ret(this->Size());
		for (int c = 0; c < this->comp.size(); c++)
			ret.comp[c] = this->comp[c]*d;

		return ret;
	}

	template <typename T>
	Vector <T>
	Vector<T>::operator+ (const Vector &v)
	{
	    Vector <T> ret(this->Size());
		for (int c = 0; c < this->comp.size(); c++)
			ret.comp[c]=this->comp[c] + v[c];

		return ret;
	}

	template <typename T>
	Vector <T> &
	Vector<T>::operator+= (const T &v)
	{
		for (int c = 0; c < this->comp.size(); c++)
			this->comp[c] += v;

		return *this;
	}



	template <typename T>
	const std::string Vector<T>::outstr() const
	{

		std::string s;
		std::ostringstream strs;
		//strs << "elems"<<this->Size()<<endl;
		for (int i = 0; i < this->Size(); i++)
			strs << std::setprecision(3)<<this->comp[i] << " ";

		s+= strs.str();
		s += " \n";

		return s;
	}

}//FluxSol

#endif
