#ifndef SistCoord_H
#define SistCoord_H

#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>

using namespace std;

class SistCoord{

	public:
	SistCoord(int i, int tipo);
	const int VerTipo();
	const int VerId();

	private:
	int id;
	int Tipo;	//0 es rect, 1 es cil y 3 es esferico

};

#endif
