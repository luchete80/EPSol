#ifndef _SINGLE_INPUT_FILE_H_
#define _SINGLE_INPUT_FILE_H_

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdlib.h> //atoi

#include "FEGrid.h"
#include "Material.h"

#include "Matrix.h" //Data Blocks

using namespace FluxSol;

template <typename number>
const std::vector <number> ReadData(string data, int *pos, string separator);
//Reads comma separated data until end of line
const bool IsValidNumber(string data, int posi, string separator);
void strip_white_spcs(std::string str);

//TO MODIFY
//TEMPLATIZE TYPE AND INHERITED LOAD AND BC FROM ONE ParENT CONDITION
//TO MODIFY, SWITCH LOADS AND BOUNDARIES INPUT TO FIELD DIM
class BoundaryFix
{
	std::vector<bool> fix;
	std::vector<double> value;  //Fixity value: displacements, temperatures, etc
	int nodeid;

	public:
		BoundaryFix(){};
		BoundaryFix(const int &id, bool x, bool y)
		{
			nodeid=id;
			fix.push_back(x);fix.push_back(y);
		};

		void BoundaryValue(double x, double y)
		{
			value.push_back(x);value.push_back(y);
		}
		const int & NodeId()const{return nodeid;}
		const std::vector<bool> & Fix ()const{return fix;}
		const std::string Show() const
		{
			string cad;
			std::ostringstream strs;
			strs<<nodeid;
			cad+=strs.str();
			cad+="\n";
			for (int f=0;f<fix.size();f++)
			{
				strs.clear();strs<<fix[f]<<";";
				cad+=strs.str();
				cad+="\n";
			}
			const string ret(cad);
			return ret;
		}

		const std::vector<double> & Value ()const{return value;}


};

class LoadC
{
	std::vector<double> load;
	int nodeid;

	public:
		LoadC(){};
		LoadC(const int &id, const double &x, const double &y)
		{
			nodeid=id;
			load.push_back(x);load.push_back(y);
		}

		const int & NodeId()const{return nodeid;}
		const std::vector<double> & Load ()const{return load;}

};

template <int dim>
class SingleInputFile
{

protected:
    string rawData;
    fstream file;

    //const int dim;
	//Temporary
	std::vector < FluxSol::Element<dim> *> ve;	//TO MODIFY DIMENSION!!!

	std::vector < FluxSol::Node > vn;

	FluxSol::FeGrid<dim> fegrid;

	double E,nu;

	FluxSol::Material * mat;	//TO MODIFY Convert to vector of pointers

	std::vector<BoundaryFix> bfix;
	std::vector<LoadC> loads;

	int bfixnumber;
	int loadnumber;

	public:
	//Constructors
	SingleInputFile(){}
	SingleInputFile(const string &name);

	void strip_white_spaces();
	void strip_endlines();

	void ShowData();
	void ReadNodes();
	void ReadElements();
	void ReadMaterial();
	void ReadBoundary();
	void ReadLoads();
	void ReadSections();

	template <int d>
	const Element <d> Function();
	const int & BFixNumber()const{return bfixnumber;}
	const BoundaryFix & BFix(const int &i)const{return this->bfix[i];}

	const int & BLoadNumber()const{return loadnumber;}
	const LoadC & BLoad(const int &i)const{return this->loads[i];}

	bool extract_in_between(string &data, string begin, string end, string &result,bool check_char_before, string acceptList);
	const FluxSol::FeGrid<dim> & Grid() const {return this->fegrid;}

	const std::vector<FluxSol::Element<dim> *> & Ve()const{return this->ve;}
	const std::vector<FluxSol::Node > & Vn()const{return this->vn;}


	const Material & Mat()const {return *this->mat;}	//TO MODIFY, MATERIAL VECTOR
};


//Reads comma separated data until end of line
template <typename number>
const std::vector <number> ReadData(string data, int *pos, string separator)
{

	vector<number> ret;

	int found=0;
	//TO MODIFY
	int foundend=data.find("\n",*pos);		//End of nodes
	//if (foundend==string::npos)

	string cut;
	string hola;
	vector <number> vdata;

	int i=0;
	int add;

	//cout << "foundend"<<endl;

	if (*pos<foundend)
	{
		//cout <<"First string char to be analyzed: " <<data[*pos] <<endl;

		bool end=false;
		while(!end)
		{
			int posint=data.find(",",*pos);
			int posend=data.find("\n",*pos);
			add=0;
			if (posend<posint || posint==string::npos)
			{
				posint=posend;
				end=true;
				add=1;
			}
			//cout<<"Next comma pos: "<<posint<<endl;
			cut = data.substr(*pos,posint-*pos+1-add);

			//number num=atoi(cut.c_str());
			char  *pend;
			number num = strtod(cut.c_str(),&pend);
			ret.push_back(num);
			//cout <<"cut string: "<<cut<<endl;
			//cout <<"value read: " << num <<endl;
			*pos= posint +1;
			i++;
		}
	}

	const vector<number> cret(ret);
	return cret;

}



//Reads comma separated data until end of line
template <typename number>
const std::vector < std::vector < number > > ReadDataBlock(string data, int *pos, string separator)
{

	vector < vector<number> > ret;

	int found=0;
	//TO MODIFY
	int foundend=data.find("\n",*pos);		//End of nodes
	//if (foundend==string::npos)

	string cut;
	string hola;
	vector <number> vdata;

	int i=0;
	int add;

	//cout << "foundend"<<endl;


	bool end=false;
	cout << "string size" <<data.length() << endl;
	bool valid=false;
	while(!end)
	{

		if (IsValidNumber(data,pos,separator))
		{
			//cout << "reading data at pos " << *pos <<endl;
			//cout << "is valid number"<<endl;
			vector <number> line=ReadData<number>(data,pos,separator);
			//cout << "pushing back line"<<endl;
			ret.push_back(line);

			//cout << "pos " <<*pos <<endl;
			//cout << "Char: "<< data.c_str()[*pos] << endl;
			if (*pos>=data.length())	{end=true; cout << "Ended"<<endl;}
			//pos is now located in first char from next line
			valid=true;
		}
		else	{cout <<"Not valid"<<endl;end=true;}
	}
	//cout << "End of reading " <<endl;
	number n=0;

	const vector < vector <number> > cret(ret);

	if (!valid)
		cret=vector < vector <number> > (n);

	else
	cret=ret;
	//cout << "Vector created"<<endl;
	return cret;

}

template <typename number>
const Matrix< number > ReadDataBlockM(string data, int *pos, string separator)
{
	vector <vector <number> > v= ReadDataBlock<number>(data,pos,separator);

	//cout << "created"<<endl;
	int r=v.size();
	int c=0;
	if (r>0)	c=v[0].size();

	//cout << "Creating matrix"<<endl;
	Matrix <number> m(r,c);

	//cout << "rows cols" << r << " " <<c<<endl;
	for (int i=0;i<r;i++)
		for (int j=0;j<v[i].size();j++)
			m[i][j]=v[i][j];

	Matrix<number> cmr(m);
	//cout << "Block" <<cmr.outstr();
	return cmr;

}

//Finds last valid position of integers
template <typename number>
const int LastValidNumberPos(string data, int *pos)
{
	int ret;
	string cut;
	int add;

	bool end=false;
	while(!end)
	{

		int foundend=data.find("\n",*pos);		//End of nodes

		if (foundend==string::npos)
			end=true;							//End of file found
		else
		{
			int posint=data.find(",",*pos);
			int posend=data.find("\n",*pos);

			if (posend<posint || posint==string::npos)
			{
				posint=posend;
				end=true;

			}
			//cout<<"Next comma pos: "<<posint<<endl;
			cut = data.substr(*pos,posint-*pos+1-add);

			char *pend;
			number num = strtod(cut.c_str(),&pend);

			if (num==0.0)	end=true;
			else 			*pos= posint +1;
		}
	}

	cout << "Next valid position" <<endl;
	return ret;

}

#endif
