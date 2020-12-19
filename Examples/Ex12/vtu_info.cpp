	file << "<?xml version=\"1.0\"?>" << endl;
	file << "<VTKFile type=\"UnstructuredGrid\">" << endl;
	file << "<UnstructuredGrid>" << endl;
	file << "<Piece NumberOfPoints=\"" << grid.Num_Verts() << "\" NumberOfCells=\"" << grid.Num_Cells() << "\">" << endl;
	file << "<Points>" << endl;
	file << "<DataArray NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\" >" << endl;
	for (int n=0;n<grid.Num_Verts();++n) {
		for (int i=0; i<3; ++i) file<< setw(16) << setprecision(8) << scientific << grid.Vertex(n).Comp()[i] << endl;
	}
	file << "</DataArray>" << endl;
	file << "</Points>" << endl;
	file << "<Cells>" << endl;
	
	file << "<DataArray Name=\"connectivity\" type=\"Int32\" format=\"ascii\" >" << endl;
	for (int c=0;c<grid.Num_Cells();++c) {
		for (int n=0;n<grid.Cell(c).Num_Vertex();++n) {
			file << grid.Cell(c).Vert(n) << "\t";
		}
		file << endl;
	}
	
	file << "</DataArray>" << endl;
	file << "<DataArray Name=\"offsets\" type=\"Int32\" format=\"ascii\" >" << endl;
	int offset=0;
	for (int c=0;c<grid.Num_Cells();++c) {
		offset+=grid.Cell(c).Num_Vertex();
		file << offset << endl;
	}
	file << "</DataArray>" << endl;
	
	file << "<DataArray Name=\"types\" type=\"UInt8\" format=\"ascii\" >" << endl;
	for (int c=0;c<grid.Num_Cells();++c) {
		if (grid.Cell(c).Num_Vertex()==4) file << "10" << endl; // Tetra
		if (grid.Cell(c).Num_Vertex()==8) file << "12" << endl; // Hexa
		if (grid.Cell(c).Num_Vertex()==6) file << "13" << endl; // Prism
		if (grid.Cell(c).Num_Vertex()==5) file << "14" << endl; // Pyramid (Wedge)
	}
	file << endl;
	file << "</DataArray>" << endl;;
	
	file << "</Cells>" << endl;
	
	file << "<CellData Scalars=\"scalars\" format=\"ascii\">" << endl;
	
	//Begin data field output
	file << "<DataArray Name=\"";

	file << "Var";
	file << "\" type=\"Float32\" format=\"ascii\" >" << endl;

	//If scalars
	for (int n=0;n<grid.Num_Cells();n++)
		file << field.Val(n).Comp()[0] << endl;

	file << "</DataArray>" << endl;

	// End of data output
	file << "</CellData>" << endl;
	
	file << "</Piece>" << endl;
	file << "</UnstructuredGrid>" << endl;
	file << "</VTKFile>" << endl;
	file.close();