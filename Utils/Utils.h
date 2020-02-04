#ifndef VARIOS_H
#define VARIOS_H

//ACA NO TENGO QUE INCLUIR A MODELO PORQUE LO USA FVGRID; TENGO QUE MANTENERLO APARTE A NASTRAN DE FVGRID PARA QUE NO SE LLAMEN SOLAS
//

#include <sstream>
#include <limits>
#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <stdio.h>
//#include <initializer_list>

using namespace std;

string int2str(int number);

int Leer_Campo(const std::string cad,const int campo);

template <typename T>
void listar_vec(const std::vector <T>, ofstream f);
//Estas funciones solo buscan la primera posicion del vector (o la unica)
int & SearchVal(const int &i,std::vector <int> &v);
bool  FindAllVals(std::vector <int> &ind, std::vector <int> &v);	//Devuelve true si encuentra todos
//Deberia ser inline
int & NumVecElemFound(std::vector <int> &ind, std::vector <int> &v);

std::vector <int> & CommonValues(std::vector <int> &ind, std::vector <int> &v);
bool ValueFound(const int &val, std::vector<int> vec);

//POnerlo como inline
template <typename T,unsigned S>
unsigned arraysize(/*const*/ T (&v)[S]) { return S; }

//template <typename T>
//class MyVector
//{
//private:
//    vector<T> m_Array;
//
//public:
//    MyVector() { }
//
//    MyVector(const initializer_list<T>& il)
//    {
//        // Manually populate the elements of the array from initializer_list x
//        for (auto x: il) // use range-based for statement to iterate over the elements of the initializer list
//            m_Array.push_back(x); // push them into the array manually
//    }
//};
//
//template <typename T>
//class MyList
//{
//private:
//    list<T> m_Array;
//
//public:
//    MyList() { }
//
//    MyList(const initializer_list<T>& il)
//    {
//        // Manually populate the elements of the array from initializer_list x
//        for (auto x: il) // use range-based for statement to iterate over the elements of the initializer list
//            m_Array.push_back(x); // push them into the array manually
//    }
//};

//const int & max(const int &a, const int &b)
//{
//	if (a>b)return a;
//	else return b;
//}



#endif
