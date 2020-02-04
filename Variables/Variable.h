#ifndef VARIABLE_H_
#define VARIABLE_H_


#include <vector>

#include "../Type/Vec3D.h"
#include "../Type/Scalar.h"
//#include <boost/variant.hpp>

namespace FluxSol
{


	//Puedo probar esto
	//O tambien las variables pueden ser derivadas
	// Por ejemplo de vectorial y escalar, pero por ejemplo, la variable puede ser una clase
	//eso lo hace mucho mas general
	//SCALAR, ETC
	template <class T>
	class Variable{

		int Id;     //Id de la variable
		T value;

	public:
		//Constructores
		Variable(){};
		Variable(const T &){};


		//USED AS CONSTRUCTION ONLY
		//RELATED TO TYPE
		//virtual const int dim(){ return 0; }

		T & Valor(){ return value; };


	};


	class TempVar
		:public Variable<Scalar>
	{


		public:


	};


	class PressVar
		:public Variable<Scalar>
	{

	};

	class DispVar
		:public Variable<Vec3D>
	{



	};

	class RotVar
	:public Variable<Vec3D>
	{



	};




	//ELEMENTS HAVE THIS
	class VarScheme{

		//instead of smart base pointers

		//typedef boost::variant <PressVar, TempVar> var_variants;

		//IR DEFINING EVERYONE
		PressVar pv;

		//std::vector<Variable <T>*> v;
		//std::vector<var_variants>;
		//
		std::vector <int> var_dof;

		public:


	};


	class DispScheme :
		public VarScheme
	{



	};

}

#endif // VARIABLE_H_INCLUDED
