#ifndef CADENAS_H
#define CADENAS_H

#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <iomanip>
#include <locale>
#include <sstream>
#include <math.h>

using namespace std;

string cadena_alineada(int tam,int val);
const double stringadouble(string Text);
const string dobleastring_alineado (const double &val);	//Tipo 0 es un double 1 es un int
const string intastring_alineado (const int &val);
string Leer_Hasta_Caracter(const string cad, char ch);

#endif