	cout << "METIS..." <<endl;
	//PETSC_Solver<double,2> solver(dofhandler.NumDoF());
	//TO MODIFY
	//Can be an integer or a non constant vector
	//solver.PreAllocateRows(dofhandler.Adj_DoF_Number());
	
	CSRGrid <2> csr(grid);
	//cout << csr.outstr();
	
	//With own CSR
	idx_t *perm, *iperm;
	
	int err;
	//err=m2gmetis(argc, args	);
	
	
	//With METIS conversion
	
	//int METIS MeshToNodal(idx t *ne, idx t *nn, idx t *eptr, idx t *eind, idx t *numflag,
	//idx t **xadj, idx t **adjncy)

	idx_t *xadj,*adjncy;
	
	std::vector <idx_t> eptr, eind;
	
	int eindsize=0;
	for (int e=0;e<grid.NumElem();e++)
		eindsize+=grid.Elem(e).NumNodes();
	
	idx_t *eind2, *eptr2;
	
	cout << "eindsize" << eindsize << endl;
	
	//METIS_Free();
	
	eptr2 = new idx_t[eindsize];
	eind2 = new idx_t[grid.NumElem()+1];
	
	cout << "allocating vectors..."<<endl;
	int pos=0;
	eptr.push_back(0);
	int poseind=0;
	for (int e=0;e<grid.NumElem();e++)
	{
		Element <2> el(grid.Elem(e));
		for (int ne=0;ne<el.NumNodes();ne++)
		{
			eind.push_back(el.NodePos(ne));
			eind2[poseind]=el.NodePos(ne);
			poseind++;
		}
		pos+=el.NumNodes();
		eptr.push_back(pos);
		eptr2[e+1]=pos;			//eptr stores the begin of each element to eind
	}
	
	cout <<"eptr"<<endl;
	for (int n=0;n<eptr.size();n++)
		cout << eptr2[n] <<endl;

	cout <<"eind"<<endl;
	for (int n=0;n<eind.size();n++)
		cout << eind2[n] <<endl;
		
	cout << "Converting Mesh to Nodal... "<<endl;
	//int METIS MeshToNodal(idx t *ne, idx t *nn, idx t *eptr, idx t *eind, idx t *numflag,
	//idx t **xadj, idx t **adjncy)	
	
	idx_t numflag=0;
	//METIS_MeshToNodal(&grid.NumElem(), &grid.NumNodes(), &eptr[0], &eind[0],&numflag,&xadj,&adjncy);

	idx_t *nv=&grid.NumNodes();
	
	idx_t *ne=&grid.NumElem();

	cout << "Num elem " << *ne;
	cout << "Num nodes " << *nv;

	
	int error;
	
	for (int n=0;n<eind.size();n++) eind2[n]=eind2[n]+1;
	error = METIS_MeshToNodal(ne, nv, eptr2, eind2,&numflag,&xadj,&adjncy);
	
	//TO c++ style cast
	//METIS_MeshToNodal(&nE, &nN, const_cast<idxtype*>(&IEN[0]), &etype,
	//&pnumflag, &nxadj[0], &nadjncy[0]);
	
	//int METIS MeshToNodal(idx t *ne, idx t *nn, idx t *eptr, idx t *eind, idx t *numflag,
	//idx t **xadj, idx t **adjncy)
	
	cout << "error : " <<error<<endl;
	
	cout << "Viewing graph..." <<endl;
	
	for (int n=0;n<grid.NumNodes()+1;n++)
		cout <<"xadj" << n << " : " << xadj[n] <<endl;
		
	for (int n=0;n<24;n++)
		cout <<"adjncy" << n << " : " << adjncy[n] <<endl;
		
		
	
	
	//int METIS NodeND(idx t *nvtxs, idx t *xadj, idx t *adjncy, idx t *vwgt, idx t *options,
	//idx t *perm, idx t *iperm)
	//METIS_NodeND(&grid.NumNodes(), xadj, adjncy,NULL,NULL, perm, iperm);
				  
				  
	cout << "With own csr ..." <<endl;
	
	idx_t *xadjp=csr.XAdj();
	idx_t *adjncyp=csr.Adjncy();
	
	for (int x=0;x<	10;x++)
		cout << "xadjp " << x << ": " << xadjp[x] <<endl;
		
	for (int x=0;x<	24;x++)
		cout << "adjncy " << x << ": " << adjncyp[x] <<endl;
	
	
	idx_t options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);
	
	//
	cout << "Numbering: " <<options[METIS_OPTION_NUMBERING]<<endl;
	
	//for (int x=0;x<	10;x++)	xadjp[x]=xadjp[x];
		
	//for (int x=0;x<	24;x++) adjncyp[x]=adjncyp[x]+1;
	
	perm=new idx_t[*nv];
	iperm=new idx_t[*nv];
	
	cout << "Node ND..." <<endl;
	error=METIS_NodeND(nv, xadjp, adjncyp,NULL,options, perm, iperm);

	for (int x=0;x<	*nv;x++) cout << "perm " << x << ": " << perm [x] <<endl;
		
	for (int x=0;x<	*nv;x++) cout << "iperm" << x << ": " << iperm[x] <<endl;
	
	cout << "error : " <<error<<endl;
	
	//METIS_Free();
	delete xadj;
	delete adjncy;
	