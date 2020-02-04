#include "./Model/Model.h"

using namespace FluxSol;

Fv_CC_Model::Fv_CC_Model(const Modelo &nasmod)
{

	//Asumo que el modelo ya lo lei
	//Primero chequeo las dimensiones de la malla
	//recorro los elementos y si no hay elementos 3D la malla es bidimensional
	//recordar que elementos CQUAD puede haber en un modelo 3D;
	int e=0;
	bidim_model=true;
	//Tengo que igualar el miembro
	//nasmodel=nasmod;
	//while (e<this->numelem)
	//{
	//	if (this->Elementos[e].Tipo()=="CTETRA" || this->Elementos[e].Tipo()=="CHEXA")
	//	{			bidim_model=false;
	//		break;
	//	}
	//	e++;
	//}

	//vector <_Vertex> vv=this->Convert_VertexVector();
	//
	//Boundary bound;
	//vector <Cell_CC> vcell;
	//this->Extract_Cells_and_Boundary(vcell,bound);


	////El constructor se encarga de crear los nodos, faces y otras cuestiones que exceden a NASTRAN
	//Fv_CC_Grid fvmesh(vv,vcell,bound);


}


void Fv_CC_Model::Extract_Cells_and_BoundaryFromNastran()
{
	vector <int> patchnum;	//Numero de la propiedad
	int nump=0;
	int numboundcells=0,numboundfaces=0,numintcells=0;
	vector <Cell_CC> vc;
	vector <Patch> vp; 
	

	//Recorro todos los elementos
	//Aca tambien estan las condiciones de borde
	//for (int e=0;e<NasFile().;e++)
	//{
	//	if (bidim_model)
	//	{
	//		//Discrimino los CQUAD4 y CQUADR
	//		if (this->Elementos[e].Tipo()=="CQUAD4" || this->Elementos[e].Tipo()=="CQUADR" ||
	//			this->Elementos[e].Tipo()=="CTRIA3" || this->Elementos[e].Tipo()=="CTRIAR" )
	//		//Es un elemento del interior
	//		{
	//			//Este lo paso al vector
	//			Cell_CC cell(this->Elementos[e].VerId(),this->Elementos[e].Conect());
	//			vc.push_back(cell);
	//			numintcells++;
	//		}
	//		else //Elemento de Boundary
	//		{
	//			int currpatch;

	//			//Si las propiedades cambian son distintos patches
	//			bool enc=false;	//Una propiedad existente
	//			int p=0;
	//			
	//			while (!enc && p<patchnum.size())
	//			{
	//				if(patchnum[p]==this->Elementos[e].VerPid())
	//				{
	//					enc=true;
	//					currpatch=p;
	//				}
	//				p++;
	//			}

	//			if (!enc ||nump==0)//Si ese patch no esta	
	//			{
	//				//patchnum.push_back(this->Elementos[e].VerPid());
	//				currpatch=patchnum.size()-1;
	//				nump++;
	//				//Inserto el patch en el vector 
	//				//Patch patch(this->Elementos[e].VerPid());
	//				//vp.push_back(patch);
	//				//bound.AddPatch(patch);
	//			}
	//			
	//			//Agrego el face al patch que corresponda
	//			vp[currpatch].AddFace(numboundfaces);			
	//			numboundfaces++;
	//		}

	//	}
	//	else	//Modelo tridimensional
	//	{}
	//}//Fin del numero de elementos

	//Incluyo el numero de patches al modelo

}

