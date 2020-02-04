#include "SingleInputFile.h"

using namespace std;
using namespace FluxSol;
//using namespace FluxSol;
//Open Filename
SingleInputFile::SingleInputFile(const string &stringname)
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
	ReadLoads();
	//Reading Elements
	this->fegrid=FluxSol::FeGrid<2>(this->vn,this->ve);

	//Reading Problem type

}

void SingleInputFile::ReadNodes()
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
		cout <<"End of line"<< pos <<endl;
		cout <<"pos: "<<pos<<endl;
		cout <<"foundend: "<<foundend<<endl;

		while (pos<foundend)
		{
			cout <<"First string char to be analyzed: " <<rawData[pos] <<endl;
			int nodeid;
			double data[4];	//Id and nodes
			for (int i =0; i<3;i++)
			{
				int posint=rawData.find(",",pos);
				cout<<"Next comma pos: "<<posint<<endl;
				cut = rawData.substr(pos,posint-pos+1);
				data[i]=atoi(cut.c_str());
				cout <<"cut string: "<<cut<<endl;
				cout << data[i] <<endl;
				pos= posint +1;
			}

			//Last data does not have a comma at the end
			int posint=rawData.find("\n",pos);
			cut = rawData.substr(pos,posint);
			data[3]=atoi(cut.c_str());
			cout << data[i] <<endl;
			pos= posint +1;

			nodeid=data[0];
			this->vn.push_back(Node(nodeid,data[1],data[2],data[3]));

			//Add Node
		}

	}
}


void SingleInputFile::ReadElements()
{

	int found=0;
	string nodstring("*Element");
	//found=rawData.find_first_of(nodtring);
	found=rawData.find(nodstring);

	cout << "Reading Elements"<<endl;
	//TO MODIFY
	int foundend=rawData.find("*Material",found);		//End of nodes

	int pos;
	pos = found;
	string cut;
	string hola;
	if (found!=string::npos)	//If found *Node
	{
		pos=rawData.find("\n",found);
		pos+=1;
		cout <<"End of line"<< pos <<endl;
		cout <<"pos: "<<pos<<endl;
		cout <<"foundend: "<<foundend<<endl;

		while (pos<foundend)
		{
			cout <<"First string char to be analyzed: " <<rawData[pos] <<endl;
			int nodeid;
			std::vector<int> vdata;	//Id and nodes connectivity
			int data[5];
			vdata.assign(4,0);
			for (int i =0; i<4;i++)
			{
				int posint=rawData.find(",",pos);
				cout<<"Next comma pos: "<<posint<<endl;
				cut = rawData.substr(pos,posint-pos+1);
				//TO MODIFY, THIS MUST SEARCH ALL NODES TILL FIND POSITION
				data[i]=atoi(cut.c_str())-1;
				cout <<"cut string: "<<cut<<endl;
				cout << data[i] <<endl;
				pos= posint +1;
			}

			//Last data does not have a comma at the end
			int posint=rawData.find("\n",pos);
			cut = rawData.substr(pos,posint);
			data[4]=atoi(cut.c_str())-1;
			cout << data[i] <<endl;
			pos= posint +1;

			for (int i=0;i<4;i++)	vdata[i]=data[i+1];
			//cout <<"Creating Element"<<endl;
			QuadLinearElement elem;
			this->ve.push_back(Element<2>(vdata));

			vdata.clear();
			//Add Node
		}

	}
}
void SingleInputFile::ShowData()
{
	cout << this->rawData << endl;
}

void SingleInputFile::strip_white_spaces() {

	string whitespaces(" ");
	int found=0;
	while (found!=string::npos) {
		found=rawData.find_first_of(whitespaces);
		if (found!=string::npos) rawData.erase(found,1);
	}
}

void SingleInputFile::ReadMaterial()
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
		pos=rawData.find("\n",found);
		cout <<"Pos"<<pos;
		pos+=1;
		cout <<"End of line"<< pos <<endl;
		cout <<"pos: "<<pos<<endl;
		cout <<"foundend: "<<foundend<<endl;
		vector <double> vdata= ReadData<double>(rawData,&pos,",");

		cout<<vdata[0]<<";"<<vdata[1]<<endl;
		this->E=vdata[0];
		this->nu=vdata[1];
	}

}

void SingleInputFile::ReadBoundary()
{
	int found=0;
	string nodstring("*Boundary");
	found=rawData.find(nodstring);

	cout << "Reading Boundary: "<<endl;
	//TO MODIFY
	int foundend=rawData.find("*Loads",found);		//End of nodes

	int pos;
	pos = found;
	string cut;
	string hola;
	vector <int> vdata;

	if (found!=string::npos)	//If found *Node
	{
		pos=rawData.find("\n",found);
		pos+=1;
		cout <<"pos: "<<pos<<endl;
		cout <<"foundend: "<<foundend<<endl;
		while (pos<foundend)
		{
			vdata= ReadData<int>(rawData,&pos,",");
			cout << "Readed Size: "<<vdata.size()<<endl;
			cout << "End pos" <<pos<<endl;
			for (int i=0;i<vdata.size();i++)
			cout <<"Readed: "<<vdata[i]<<endl;;



			//TO MODIFY, SEARCH ID
			int idnode=vdata[0]-1;
			bool dir[2];dir[0]=dir[1]=false;
			//Saving Boundary Conditions
			for (int bc=1;bc<vdata.size();bc++)
			{
					dir[vdata[bc]-1]=true;
			}

			bfix.push_back(BoundaryFix(idnode,dir[0],dir[1]));
			cout << "Node: " <<idnode << "Dir 1 " << dir [0] << "Dir 2 " << dir [1]<<endl;
		}

	}
}

void SingleInputFile::ReadLoads()
{
	int found=0;
	string nodstring("*Loads");
	found=rawData.find(nodstring);

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
	cout <<"pos: "<<pos<<endl;
	cout <<"foundend: "<<foundend<<endl;
	while (pos<foundend)
	{
		vdata= ReadData<double>(rawData,&pos,",");
		cout << "Readed Size: "<<vdata.size()<<endl;
		cout << "End pos" <<pos<<endl;
		for (int i=0;i<vdata.size();i++)
		cout <<"Readed: "<<vdata[i]<<endl;;



	cout << "Pushing Back Load " <<endl;
	//Supposing one direction per line
	//TO MODIFY
	int idnode=vdata[0]-1;	//cast to int
	double dir[2];dir[0]=dir[1]=0.;
	//Saving Boundary Conditions
	dir[(int)(vdata[1])-1]=vdata[2];

	loads.push_back(LoadC(idnode,dir[0],dir[1]));

	cout << "Node: " <<idnode << "Dir 1 " << dir [0] << "Dir 2 " << dir [1]<<endl;
	}


}

void SingleInputFile::strip_endlines() {

	string whitespaces(" \t\f\v\n\r");
	int found=0;
	while (found!=string::npos) {
		found=rawData.find_first_of(whitespaces);
		if (found!=string::npos) rawData.erase(found,1);
	}
}

bool SingleInputFile::extract_in_between(string &data, string begin, string end, string &result,bool check_char_before, string acceptList) {
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
