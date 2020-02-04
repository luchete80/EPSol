#include "SingleInputFile.h"

using namespace std;
using namespace FluxSol;
//using namespace FluxSol;
//Open Filename

void strip_white_spcs(std::string str) {

	string whitespaces(" ");
	int found=0;
	while (found!=string::npos) {
		found=str.find_first_of(whitespaces);
		if (found!=string::npos) str.erase(found,1);
	}
}

template <int dim>
SingleInputFile<dim>::SingleInputFile(const string &stringname)
{
    this->file.open(stringname.c_str());
    string line;
	this->rawData="";
	while(getline(file, line)) this->rawData += line + "\n";
	file.close();

	//Material
	double nu,E;

	//
	this->strip_white_spaces();

	//Reading Nodes
	ReadNodes();
	ReadElements();
	ReadMaterial();
	ReadBoundary();
	cout << "Reading Loads"<<endl;
	ReadLoads();
	//Reading Elements
	//Taking reference from element
	cout << "Creating Mesh from data.." <<endl;
	//this->fegrid=FluxSol::FeGrid<2>(this->vn,this->ve);
    std::vector < FluxSol::Element<dim> *> v;

    this->fegrid=FeGrid<dim>(this->vn,ve);



	//Reading Problem type

}


template <int dim>
void SingleInputFile<dim>::ReadNodes()
{

	int found=0;
	string nodstring("*Node");
	//found=rawData.find_first_of(nodtring);
	found=rawData.find(nodstring);
	cout << "Reading Nodes"<<endl;
	//TO MODIFY
	int foundend=rawData.find("*Element",found);		//End of nodes

	int pos;
	pos = found;
	string cut;
	string hola;
	if (found!=string::npos)	//If found *Node
	{
		pos=rawData.find("\n",found);
		pos+=1;
		//cout <<"End of line"<< pos <<endl;
		//cout <<"pos: "<<pos<<endl;
		//cout <<"foundend: "<<foundend<<endl;

		while (pos<foundend)
		{
			//cout <<"First string char to be analyzed: " <<rawData[pos] <<endl;
			int nodeid;
			double data[4];	//Id and nodes
			for (int i =0; i<3;i++)
			{
				int posint=rawData.find(",",pos);
				//cout<<"Next comma pos: "<<posint<<endl;
				cut = rawData.substr(pos,posint-pos+1);
				char *pend;
				data[i]=strtod(cut.c_str(),&pend);
				//cout <<"cut string: "<<cut<<endl;
				//cout << data[i] <<endl;
				pos= posint +1;
			}

			//Last data does not have a comma at the end
			int posint=rawData.find("\n",pos);
			cut = rawData.substr(pos,posint);
			data[3]=atoi(cut.c_str());
			//cout << "Node Coordinates"<<endl; for (int i=1;i<4;i++) cout << data[i] <<endl;
			pos= posint +1;

			nodeid=data[0];
			this->vn.push_back(Node(nodeid,data[1],data[2],data[3]));

			//Add Node
		}

	}
	else
	{
		cout << "ERROR: No nodes defined." <<endl;
	}
}

template <int dim>
void SingleInputFile<dim>::ReadElements()
{

	int found=0;
	string nodstring("*Element");
	//found=rawData.find_first_of(nodtring);
	found=rawData.find(nodstring);
	cout <<"Element position " << found<<endl;
	cout << "Reading Elements"<<endl;

	//TO MODIFY
	int foundend=rawData.find("*Material",found);		//End of nodes, TO MODIFY


	int pos;
	pos = found;
	string cut;
	string hola;
	if (found!=string::npos)	//If found *Node
	{
		//FIND ELEMENT TYPE
		int poscomma=rawData.find(",",found);
		pos=rawData.find("\n",found);

		cut = rawData.substr(poscomma+1,pos-2-poscomma+1);
		strip_white_spcs(cut);
		//cout << "Element Type " << "-" <<cut <<"-"<< endl;

		bool eltypeok=true;

		int elementnodes;
		int elementorder;
		string eltype=cut;
		//TO MOdify, NEXT EACH CHARACTER HAS A MODE
		if (eltype=="Q4PSBF")	//Full Integration Quad 4 node Plain Stress Base element, full Integration
		{
			elementnodes=4;
			elementorder=1;
		}
		else if (eltype=="Q4PEBF")	//Full Integration Quad 4 node Plain Strain Base element
		{
			elementnodes=4;
			elementorder=1;
		}
		else if (eltype=="Q4TPLF" || eltype=="Q4TPLR") //Full Integration Quad 4 node Thermal  element
		{
			elementnodes=4;
			elementorder=1;
		}
        else if (eltype=="H8TPLF" || eltype=="H8TPLR") //Full Integration Quad 4 node Thermal  element
		{
			elementnodes=8;
			elementorder=1; //Field dim
		}
		else
		{
			eltypeok=false;
		}
		if (eltypeok)
		{
			//End reading element type
			pos+=1;
			//cout <<"End of line"<< pos <<endl;
			//cout <<"pos: "<<pos<<endl;
			//cout <<"foundend: "<<foundend<<endl;

			while (pos<foundend)
			{
				//cout <<"First string char to be analyzed: " <<rawData[pos] <<endl;
				int nodeid;
				std::vector<int> vdata;	//Id and nodes connectivity
				int data[elementnodes+1];	//Element id and element nodes
				vdata.assign(elementnodes,0);
				for (int i =0; i<elementnodes;i++)//TO MODIFY, element node index
				{
					int posint=rawData.find(",",pos);
					//cout<<"Next comma pos: "<<posint<<endl;
					cut = rawData.substr(pos,posint-pos+1);
					//TO MODIFY, THIS MUST SEARCH ALL NODES TILL FIND POSITION
					data[i]=atoi(cut.c_str())-1;
					//cout <<"cut string: "<<cut<<endl;
					//cout << data[i] <<endl;
					pos= posint +1;
				}

				//Last data does not have a comma at the end
				int posint=rawData.find("\n",pos);
				cut = rawData.substr(pos,posint);
				data[elementnodes]=atoi(cut.c_str())-1;
				//cout << data[i] <<endl;
				pos= posint +1;

				for (int i=0;i<elementnodes;i++)	vdata[i]=data[i+1];

				//TO MODIFY, insert pointer
				//TO MAKE A GENERAL VIRTUAL DATA METHOD
				if (eltype=="Q4PSBF")	//Full Integration Quad 4 node Plain Stress Base element, full Integration
				{
					Element <dim> *ep=new Q4PS<dim>(vdata);
                    this->ve.push_back(ep);
					//cout << "Pushing plain strain elements" <<endl;
				}
				else if (eltype=="Q4PEBF")	//Full Integration Quad 4 node Plain Strain Base element
				{

				}
				else if (eltype=="Q4TPLF") //Full Integration Quad 4 node Thermal  element
				{
					Element <dim> *ep=new Q4TH<dim>(vdata,2);
					this->ve.push_back(ep);
					cout << "Pushing thermal elements" <<endl;
				}
				else if (eltype=="Q4TPLR") //Full Integration Quad 4 node Thermal  element
				{
					Element <dim> *ep=new Q4TH<dim>(vdata,1);
					this->ve.push_back(ep);
					//cout << "Pushing thermal elements" <<endl;
				}
				else if (eltype=="H8TPLF") //Full Integration Quad 4 node Thermal  element
				{
					Element <dim> *ep=new H8TH<dim>(vdata,2);
					this->ve.push_back(ep);
					//cout << "Pushing thermal elements" <<endl;
				}
				else if (eltype=="H8TPLR") //Full Integration Quad 4 node Thermal  element
				{
					Element <dim> *ep=new H8TH<dim>(vdata,1);
					this->ve.push_back(ep);
					//cout << "Pushing thermal elements" <<endl;
				}
				//cout <<"Element created" <<endl;
				vdata.clear();
				//Add Node
			}
		}
		else
		{
			cout << "ERROR: Element type " << eltype <<" not recognized." <<endl;
		}//End if eltypeok
	}
	else
	{
		cout << "ERROR: No elements found." <<endl;
	}
}

template <int dim>
void SingleInputFile<dim>::ShowData()
{
	cout << this->rawData << endl;
}

template <int dim>
void SingleInputFile<dim>::strip_white_spaces() {

	string whitespaces(" ");
	int found=0;
	while (found!=string::npos) {
		found=rawData.find_first_of(whitespaces);
		if (found!=string::npos) rawData.erase(found,1);
	}
}

//TO MODIFY: PARSING COMPARES WITH SSTREAM, STRTOD, ATOF
template <int dim>
void SingleInputFile<dim>::ReadMaterial()
{
	int found=0;
	string nodstring("*Material");

	found=rawData.find(nodstring);

	cout << "Reading Material"<<endl;
	//TO MODIFY
	int foundend=rawData.find("*Boundary",found);		//End of nodes

	int pos;
	pos = found;
	string cut;
	string hola;
	size_t aver;
	if (found!=string::npos)	//If found *Node
	{
		//FIND MAETERIAL TYPE
		int poscomma=rawData.find(",",found);
		pos=rawData.find("\n",found);

		cut = rawData.substr(poscomma+1,pos-2-poscomma+1);
		strip_white_spcs(cut);
		cout << "Material Type " << "-" <<cut <<"-"<< endl;

		bool mattypeok=true;

		if (cut=="Elastic")
		{
			pos=rawData.find("\n",found);
			cout <<"Pos"<<pos;
			pos+=1;
			//cout <<"End of line"<< pos <<endl;
			//cout <<"pos: "<<pos<<endl;
			//cout <<"foundend: "<<foundend<<endl;
			vector <double> vdata= ReadData<double>(rawData,&pos,",");

			cout<<vdata[0]<<";"<<vdata[1]<<endl;
			this->E=vdata[0];	//TO MODIFY
			this->nu=vdata[1];
			LinearElasticMaterial *m=new LinearElasticMaterial(2,E, nu);
			//LinearElasticMaterial m(E, nu);
			this->mat= m;	//Direction
		}
		else
		if (cut=="Thermal")
		{
			pos=rawData.find("\n",found);
			cout <<"Pos"<<pos;
			pos+=1;
			//cout <<"End of line"<< pos <<endl;
			//cout <<"pos: "<<pos<<endl;
			//cout <<"foundend: "<<foundend<<endl;
			vector <double> vdata= ReadData<double>(rawData,&pos,",");

			cout<<vdata[0]<<";"<<vdata[1]<<endl;
			double k=vdata[0];
			LinearElasticMaterial *m=new LinearElasticMaterial(1,k); //TO MODIFY
			//LinearElasticMaterial m(E, nu);
			this->mat= m;	//Direction
		}
		else
		{
			mattypeok=false;
		}

	}

}

//Reads boundary conditions
template <int dim>
void SingleInputFile<dim>::ReadBoundary()
{
	int found=0;
	string nodstring("*Boundary");
	found=rawData.find(nodstring);

	cout << "Reading Boundary: "<<endl;
	//TO MODIFY
	int foundend=rawData.find("*Loads",found);		//End of nodes
	int posst;

	posst=rawData.find("\n",found);
	cout <<"Pos"<<posst;
	posst+=1;
	//int fouLast=LastValidNumberPos<int>(rawData,&posst);
	//vector <vector <double> >block=ReadDataBlock<double>(rawData,&posst,",");
	Matrix <double> block=ReadDataBlockM<double>(rawData,&posst,",");

	int bcnum=0;

	for (int i=0;i<block.Rows();i++)
	{
		//TO MODIFY, SEARCH INDEX
		int idnode=block[i][0]-1;
		bool dir[2];dir[0]=dir[1]=false;
		double value[2]; value[0]=value[1]=0.;

        //TO MODIFY
		//If dimension not specified, asuming all directions all these value
		if ((block.Cols())==2)
		{
		    dir[0]=dir[1]=true;
		    value[0]=value[1]=block[i][1];
		}
        else //Only fixs 1 direction
        {
            //cout << "values " << value [0] << " " <<value[1]<<endl;
            int idir=block[i][1]-1;
            dir[idir]=true;
            value[idir]=block[i][2];
        }

        BoundaryFix bf(idnode,dir[0],dir[1]);
        bf.BoundaryValue(value[0],value[1]);
		bfix.push_back(bf);
		//cout << "Node: " <<idnode << "Dir 1 " << dir [0] << "Dir 2 " << dir [1]<<endl;
		bcnum++;
	}

	this->bfixnumber=bfix.size();
	cout << "Number of Boundary conditions: " <<bcnum<<endl;
	cout << "Exiting BCS..."<<endl;
}

template <int dim>
void SingleInputFile<dim>::ReadLoads()
{
	int found=0;
	string nodstring("*Loads");
	found=rawData.find(nodstring,1);

	cout << "Reading Loads: "<<endl;
	//TO MODIFY
	//int foundend=rawData.find("*Anything",found);		//End of nodes
	int foundend=rawData.size();

	int pos;
	pos = found;
	string cut;
	string hola;
	vector <double> vdata;

	pos=rawData.find("\n",found);
	pos+=1;
	//cout <<"pos: "<<pos<<endl;
	//cout <<"foundend: "<<foundend<<endl;
	while (pos<foundend)
	{
		vdata= ReadData<double>(rawData,&pos,",");
		//cout << "Readed Size: "<<vdata.size()<<endl;
		//cout << "End pos" <<pos<<endl;
		//for (int i=0;i<vdata.size();i++)  cout <<"Readed: "<<vdata[i]<<endl;;



        //cout << "Pushing Back Load " <<endl;
        //Supposing one direction per line
        //TO MODIFY
        int idnode=vdata[0]-1;	//cast to int
        double dir[2];dir[0]=dir[1]=0.;
        //Saving Boundary Conditions
        dir[(int)(vdata[1])-1]=vdata[2];

        loads.push_back(LoadC(idnode,dir[0],dir[1]));

        //cout << "Node: " <<idnode << "Dir 1 " << dir [0] << "Dir 2 " << dir [1]<<endl;
	}

	this->loadnumber=loads.size();

}

template <int dim>
//After reading elements and materials
void SingleInputFile<dim>::ReadSections()
{


}



///////////// GENERAL ///////////////
template <int dim>
void SingleInputFile<dim>::strip_endlines() {

	string whitespaces(" \t\f\v\n\r");
	int found=0;
	while (found!=string::npos) {
		found=rawData.find_first_of(whitespaces);
		if (found!=string::npos) rawData.erase(found,1);
	}
}
template <int dim>
bool SingleInputFile<dim>::extract_in_between(string &data, string begin, string end, string &result,bool check_char_before, string acceptList) {
	int begin_pos, end_pos;
	begin_pos=0; end_pos=0;
	string pre;
	bool found=false;
	while (!found) {
		// Find the first occurance position of the beginning sequence
		begin_pos=data.find(begin,end_pos);
		if (begin_pos==string::npos) return false;
		// From where the first search left off, find the first occurance position of the ending sequence
		end_pos=data.find(end,begin_pos+begin.length());
		if (end_pos==string::npos) return false;
		// Check the character just before the beginning delimiter (if asked for)
		if (check_char_before) {
			if (begin_pos==0) {
				found=true;
			} else {
				pre=data.substr(begin_pos-1,1);
				if (pre.find_first_of(acceptList)!=string::npos) found=true;

			}
		} else {
			found=true;
		}
		if (found) {
			// Extract the what's in between
			result=data.substr(begin_pos+begin.length(),end_pos-begin_pos-begin.length());
			// Remove that chunk from the originial data
			data.replace(begin_pos,end_pos+end.length()-begin_pos,"");
		}
	}
	return true;
}

const bool IsValidNumber(string data, int posi, string separator)
{
	int *pos=posi;


	int found=0;
	//TO MODIFY
	int foundend=data.find("\n",*pos);		//End of nodes
	//if (foundend==string::npos)

	string cut;

	int i=0;
	int add;

	//cout << "foundend"<<endl;
	bool ret=false;
	if (*pos<foundend)
	{
		//cout <<"First string char to be analyzed: " <<data[*pos] <<endl;

			int posint=data.find(",",*pos);
			int posend=data.find("\n",*pos);
			add=0;
			if (posend<posint || posint==string::npos)
			{
				posint=posend;
				add=1;
			}
			//cout<<"Next comma pos: "<<posint<<endl;
			cut = data.substr(*pos,posint-*pos+1-add);

			//number num=atoi(cut.c_str());
			char  *pend;
			double num = strtod(cut.c_str(),&pend);
			if (num!=0) ret=true;
			//cout <<"cut string: "<<cut<<endl;
			//cout <<"value read: " << num <<endl;

	}
	const bool cret=ret;
	return cret;

}

#include "SingleInputFile.inst"
