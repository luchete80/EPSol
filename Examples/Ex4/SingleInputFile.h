#ifndef _SINGLE_INPUT_FILE_H_
#define _SINGLE_INPUT_FILE_H_

#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include "EPSol.h"

template <typename number>
const std::vector <number> ReadData(string data, int *pos, string separator);

//TO MODIFY
//TEMPLATIZE TYPE AND INHERITED LOAD AND BC FROM ONE ParENT CONDITION
class BoundaryFix
{
	std::vector<bool> fix;
	int nodeid;
	
	public:
		BoundaryFix(){};
		BoundaryFix(const int &id, bool x, bool y)
		{
			nodeid=id;
			fix.push_back(x);fix.push_back(y);
		};
		
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

class SingleInputFile
{

protected:
    string rawData;
    fstream file;
	
	//Temporary
	std::vector < FluxSol::Element<2> > ve;
	std::vector < FluxSol::Node > vn;

	FluxSol::FeGrid<2> fegrid;
	double E,nu;
	
	std::vector<BoundaryFix> bfix;
	std::vector<LoadC> loads;
	
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
	const int & BFixNumber()const{return this->bfix.size();}
	const BoundaryFix & BFix(const int &i)const{return this->bfix[i];}
	
	const int & BLoadNumber()const{return this->loads.size();}
	const LoadC & BLoad(const int &i)const{return this->loads[i];}
	
	bool extract_in_between(string &data, string begin, string end, string &result,bool check_char_before, string acceptList);
	const FluxSol::FeGrid<2> & Grid() const {return this->fegrid;}
	
	const std::vector<FluxSol::Element<2> > & Ve()const{return this->ve;}
	const std::vector<FluxSol::Node > & Vn()const{return this->vn;}
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
	
	cout << "foundend"<<endl;
	
	if (*pos<foundend)
	{
		cout <<"First string char to be analyzed: " <<data[*pos] <<endl;
	
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
			cout<<"Next comma pos: "<<posint<<endl;
			cut = data.substr(*pos,posint-*pos+1-add);
				
			number num=atoi(cut.c_str());
			ret.push_back(num);
			cout <<"cut string: "<<cut<<endl;
			cout <<"value read: " << num <<endl;
			*pos= posint +1;
			i++;
		}		
	}

	const vector<number> cret(ret); 
	return cret;

}

#endif
