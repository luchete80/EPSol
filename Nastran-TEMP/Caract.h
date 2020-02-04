#ifndef CARACT_H
#define CARACT_H

#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <iomanip>
#include <locale>
#include <sstream>

using namespace std;
//Grupo de caracteristicas hasta ahora la precarga y la friccion
class Caract{
	public:
		Caract();
		Caract(const double &mu,const double &po);
		const double VerMu();
		const double VerPo();
	private:
		double Mu,Po;

};


#endif
