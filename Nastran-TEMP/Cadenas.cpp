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
//Rutina para escribir un campo entero en pos, con un tamaño de caracter
//tam
#include "Cadenas.h"

const double stringadouble(string Text)
{
	double Number;
	if ( ! (istringstream(Text) >> Number) ) Number = 0.0;

	return Number;
}


//PASO UN DOBLE A UN STRING EN FORMATO NASTRAN
const string dobleastring_alineado (const double &val)
{
	char cad[100],cad2[100];
	bool fin=false;
	//Donde defino caddef??
	string cadtemp,caddef;
	caddef="";
	int decenas=(int)log10(fabs(val));
	int decenas_pos;	//decenas positivas
	if (decenas<0)	decenas_pos=0;
	else					decenas_pos=decenas;

	//cad 2 la uso para control, para no modificar este algoritmo
	sprintf(cad2,"%1.20e",val);
	string cad2s=cad2;

	int j=cad2s.size()-6;	//Retrocedo a la altura de un digito antesl signo del exponente
	while(!fin)
	{
		if (cad2s[j]!='0')
			fin=true;
		else
			j--;
	}

	int pruvvv=cad2s.find_last_of("0",1);
	//CANTIDAD DE DIGITOS SIGNIFICATIVOS, NO INCLUYEN LA COMA, OBVIAMENTE EL SIGNO SI
	int digitos_signif=j;		//Sumo el numero de antes de la coma Posicion de la coma .

	//No tiene sentido tener tanta presicion con numeros tan grandes
	//Si NO es un numero menor a 10 con muchos decimales significativos
	//Ya si es del orden de lso miles con 3 decimales ya tengo 8
	if (digitos_signif<5 || fabs(val)>1000  || fabs(val)<1.0E-03)	//ESCRIBO EN NOTACION CIENTIFICA
	{
		sprintf(cad,"%1.3e",val);

		cadtemp=cad;
		//Saco los decimales dps de la coma y
		//la cadena e-00
		int inc=0;
		//Primero busco las posiciones de la coma y de la e
		int poscoma=0;
		//La posicion de la coma dep}e de si es positivo o
		//negativo
		fin=false;
		while (!fin)
		{
			if (cadtemp[poscoma]=='.')
				fin=true;
			else
				poscoma=poscoma+1;
		} //while

		for (int i=0;i<poscoma;i++){
			char c=cadtemp[i];
			caddef+=c;}
		caddef+='.';

	//	caddef[poscoma]+=cadtemp[poscoma];
		//Si lo que está a la derecha de la coma es distinto de 0
		inc=3;  //Son 3 decimales, me fijo cual es el ultimo significativo

		//Me fijo si hay ceros dps de la coma
		int dec;
		for (dec=3;dec>0;dec--)
		{
			if (cadtemp[poscoma+inc]=='0')
				inc=inc-1;
			else
				break;

		}//Fin del for

		//inc son os decimales significativos
		for (dec=1;dec<inc+1;dec++)
			caddef+=cadtemp[poscoma+dec];



		int expsig=0;   //Cifras significativas del exponente
		for (int exp=1;exp<4;exp++)
		{
			//le sumolos 3 deciumales dps de la coma, la e y el sugno del exponente
			if (cadtemp[poscoma+3+1+1+exp]!='0')
				expsig=expsig+1;
		}

		//si hay exponente con decimales significativos
		int comienzo=3-expsig;
		if (expsig>0)
		{
			//Signo del exponente
			caddef+=cadtemp[poscoma+3+1+1];
			//le sumolos 3 deciumales dps de la coma, la e y el sugno del exponente
			for (int i=1;i<expsig+1;i++)	caddef+=cadtemp[poscoma+3+1+1+comienzo+i];
		}
	}
	else//Los digitos significativos son mayores a 5, entonces lo paso como
	{
		int decimales=digitos_signif;	//Le saco el primer numero que va antes de la coma
															//Recordando que la coma no formaba parte
		if (fabs(val)>=1.0)				//Si el modulo del numero es mayor o igual que 1
			decimales--;
		if (val<0) decimales--;
		int maxdec=8-2-decenas_pos;		//Maximos decimales si el numero es positivo
																		//Resto el primer numero y la coma
		if (val>0 && decimales>maxdec)
			decimales=maxdec;
		else if (val<0 && decimales>(maxdec-1))
			decimales=maxdec-1;
		string formato="%1.";	//para el sprintf
		char pru[100];
		sprintf(pru,"%d",decimales);
		formato+=pru;
		formato+='f';
		sprintf(pru,formato.c_str(),val);
		caddef=pru;
		//Depende de si es positivo o negativo el numero tomo mas cifras

	}//Fin de si los digitos significativos son mayores o iguales a 4
	//Ajusto a derecha la cadena
	string cadalin="        ";
	int pos=8-caddef.size();
	int dig=0;
	for (int k=0;k<caddef.size();k++)
	{
		cadalin[k+pos]=caddef[dig];
		dig=dig+1;
	}//De la copia de caracteres

	return cadalin;
}


const string intastring_alineado (const int &val)
{
	//Donde defino caddef??
	string caddef;

	char cad[100];
	sprintf(cad,"%d",val);

	caddef=cad;
	//Ajusto a derecha la cadena
	string cadalin="        ";
	int pos=8-caddef.size();
	int dig=0;
	for (int k=0;k<caddef.size();k++)
	{
		cadalin[k+pos]=caddef[dig];
		dig=dig+1;
	}//De la copia de caracteres

	return cadalin;
}

string Leer_Hasta_Caracter(const string cad, char ch)
{
	string r;

	int c=0;
	bool fin=false;
	while (c<cad.size())
	{
		if (cad[c]==ch){
			fin=true;
			break;
		}
		else
		r+=cad[c];
		c++;
	}
	return r;
}