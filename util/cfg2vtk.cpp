//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) A Meo 2016. All rights reserved.
//
//-----------------------------------------------------------------------------
//
// Program to convert vampire cfg files to vtk format (for programs such as paraview)
//
// ./cfg2vtk
//

// Standard Libraries
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

struct material_t{

	double mx,my,mz,mu_s,magm;

};

void rgb( double ireal, double &red, double &green, double &blue){

			if(ireal>0.8){
				red = 0.0;
				green = 0.0;
				blue = 1.0;
			}
			else if(ireal>=0.0){
				red = 1.0-ireal*1.2;
				green = 1.0-ireal*1.2;
				blue = 1.0;
			}
			else if(ireal>=-0.8){
				red = 1.0;
				green = 1.0+ireal*1.2;
				blue = 1.0+ireal*1.2;
			}
			else if(ireal<-0.8){
				red = 1.0;
				green = 0.0;
				blue = 0.0;
			}
			else{
				red = 1.0;
				green = 1.0;
				blue = 1.0;
			}

			if(blue<0.0) blue=0.0;
			if(red<0.0) red=0.0;
			if(green<0.0) green=0.0;

}

int main(){

	std::vector <int> mat(0);
	std::vector <int> cat(0);
	std::vector <double> spinx(0);
	std::vector <double> spiny(0);
	std::vector <double> spinz(0);
	std::vector <double> cx(0);
	std::vector <double> cy(0);
	std::vector <double> cz(0);
	std::vector <std::string> type(0);
	std::vector <std::string> filenames(0);

	// open coordinate file
	std::ifstream coord_file;
	coord_file.open("atoms-coords.cfg");

	// read in file header
	std::string dummy;

	getline(coord_file,dummy);
	//std::cout << dummy << std::endl;

	getline(coord_file,dummy);
	//std::cout << dummy << std::endl;

	getline(coord_file,dummy);
	//std::cout << dummy << std::endl;

	getline(coord_file,dummy);
	//std::cout << dummy << std::endl;

	getline(coord_file,dummy);
	//std::cout << dummy << std::endl;

	// get number of atoms
	unsigned int n_atoms;
	getline(coord_file,dummy);
	//std::cout << dummy << std::endl;
	dummy.erase (dummy.begin(), dummy.begin()+17);
	n_atoms=atoi(dummy.c_str());

	getline(coord_file,dummy);
	//std::cout << dummy << std::endl;

	// get number of subsidiary files
	unsigned int n_files;
	getline(coord_file,dummy);
	//std::cout << dummy << std::endl;
	dummy.erase (dummy.begin(), dummy.begin()+22);
	n_files=atoi(dummy.c_str());

	for(int file=0; file<n_files; file++){
		getline(coord_file,dummy);
		filenames.push_back(dummy);
		//std::cout << filenames[file] << std::endl;
	}

	getline(coord_file,dummy);
	//std::cout << dummy << std::endl;

	unsigned int n_local_atoms;
	getline(coord_file,dummy);
	//std::cout << dummy << std::endl;
	n_local_atoms=atoi(dummy.c_str());

	// resize arrays
	mat.resize(n_atoms);
	cat.resize(n_atoms);
	spinx.resize(n_atoms);
	spiny.resize(n_atoms);
	spinz.resize(n_atoms);
	cx.resize(n_atoms);
	cy.resize(n_atoms);
	cz.resize(n_atoms);
	type.resize(n_atoms);

	unsigned int counter=0;

	// finish reading master file coordinates
	for(int i=0;i<n_local_atoms;i++){
		coord_file >> mat[counter] >> cat[counter] >> cx[counter] >> cy[counter] >> cz[counter] >> type[counter];
		counter++;
	}

	// close master file
	coord_file.close();

	// now read subsidiary files
	for(int file=0; file<n_files; file++){
		std::ifstream infile;
		infile.open(filenames[file].c_str());

		// read number of atoms in this file
		getline(infile,dummy);
		n_local_atoms=atoi(dummy.c_str());
		for(int i=0;i<n_local_atoms;i++){
			infile >> mat[counter] >> cat[counter] >> cx[counter] >> cy[counter] >> cz[counter] >> type[counter];
			counter++;
		}
		// close subsidiary file
		infile.close();
	}

	// check for correct read in of coordinates
	if(counter!=n_atoms) std::cerr << "Error in reading in coordinates" << std::endl;

	// Loop over all possible spin config files
	int ios = 0;
	int spinfile_counter=0;

	while(ios==0){

		std::stringstream file_sstr;
		file_sstr << "atoms-";
		file_sstr << std::setfill('0') << std::setw(8) << spinfile_counter;
		file_sstr << ".cfg";
		std::string cfg_file = file_sstr.str();
		//const char* cfg_filec = cfg_file.c_str();

		std::ifstream spinfile;
		spinfile.open(cfg_file.c_str());

		if(spinfile.is_open()){

			std::cout << "Processing file: " << cfg_file << std::endl;

		// Read in file header
		getline(spinfile,dummy);

	//std::cout << dummy << std::endl;

	getline(spinfile,dummy);
	//std::cout << dummy << std::endl;

	getline(spinfile,dummy);
	//std::cout << dummy << std::endl;

	getline(spinfile,dummy);
	//std::cout << dummy << std::endl;

	getline(spinfile,dummy);
	//std::cout << dummy << std::endl;

	// get number of atoms
	unsigned int n_spins;
	getline(spinfile,dummy);
	//std::cout << dummy << std::endl;
	dummy.erase (dummy.begin(), dummy.begin()+17);
	n_spins=atoi(dummy.c_str());
	if(n_spins!=n_atoms) std::cerr << "Error! - mismatch between number of atoms in coordinate and spin files" << std::endl;

	getline(spinfile,dummy); // sys dimensions
	//std::cout << dummy << std::endl;
	char const field_delim = '\t';
	dummy.erase (dummy.begin(), dummy.begin()+18);
	std::istringstream ss(dummy);
	std::vector<double> val(0);
	for (std::string num; getline(ss, num, field_delim); ) {
		val.push_back(atof(num.c_str()));
		//std::cout << num << "\t" << val << std::endl;
	}
	double dim[3] = {val[0],val[1],val[2]};
	//std::cout << dim[0] << "\t" << dim[1] << "\t" << dim[2] << std::endl;

	getline(spinfile,dummy); // coord file

	getline(spinfile,dummy); // time
	//std::cout << dummy << std::endl;
	dummy.erase (dummy.begin(), dummy.begin()+5);
	double time = atof(dummy.c_str());

	getline(spinfile,dummy); // field
	//std::cout << dummy << std::endl;
	dummy.erase (dummy.begin(), dummy.begin()+6);
	double field = atof(dummy.c_str());

	getline(spinfile,dummy); // temp
	//std::cout << dummy << std::endl;
	dummy.erase (dummy.begin(), dummy.begin()+12);
	double temperature = atof(dummy.c_str());

	// magnetisation
	getline(spinfile,dummy);
	dummy.erase (dummy.begin(), dummy.begin()+14);
	std::istringstream ss2(dummy);
	val.resize(0);
	for (std::string num; getline(ss2, num, field_delim); ) {
		val.push_back(atof(num.c_str()));
		//std::cout << num << "\t" << val << std::endl;
	}
	double mx = val[0];
	double my = val[1];
	double mz = val[2];

	// get number of materials
	unsigned int n_mat;
	getline(spinfile,dummy);
	//std::cout << dummy << std::endl;
	dummy.erase (dummy.begin(), dummy.begin()+20);
	n_mat=atoi(dummy.c_str());

	std::vector <material_t> material(n_mat);
	for(int imat=0;imat<n_mat;imat++){
		getline(spinfile,dummy);
		//dummy.erase (dummy.begin(), dummy.begin()+22);
		std::istringstream ss(dummy);
		std::vector<double> val(0);
		for (std::string num; getline(ss, num, field_delim); ) {
			val.push_back(atof(num.c_str()));
		}
		material[imat].mu_s = val[0];
		material[imat].mx = val[1];
		material[imat].my = val[2];
		material[imat].mz = val[3];
		material[imat].magm = val[4];
		//std::cout << material[imat].mu_s << "\t" << material[imat].mx << "\t" << material[imat].my << "\t"  << material[imat].mz << "\t"  << material[imat].magm << std::endl;
	}

	// line
	getline(spinfile,dummy);
	//std::cout << dummy << std::endl;

	// get number of subsidiary files
	unsigned int n_files;
	getline(spinfile,dummy);
	//std::cout << dummy << std::endl;
	dummy.erase (dummy.begin(), dummy.begin()+22);
	n_files=atoi(dummy.c_str());
	filenames.resize(0);

	for(int file=0; file<n_files; file++){
		getline(spinfile,dummy);
		filenames.push_back(dummy);
		//std::cout << filenames[file] << std::endl;
	}

	getline(spinfile,dummy);
	//std::cout << dummy << std::endl;

	unsigned int n_local_atoms;
	getline(spinfile,dummy);
	//std::cout << dummy << std::endl;
	n_local_atoms=atoi(dummy.c_str());

	// Open vtk Output file
	std::stringstream vtk_file_sstr;
	vtk_file_sstr << "atoms-";
	vtk_file_sstr << std::setfill('0') << std::setw(8) << spinfile_counter;
	vtk_file_sstr << ".vtu";
	std::string vtk_file = vtk_file_sstr.str();

	std::ofstream vtkfile;
	vtkfile.open(vtk_file.c_str());

	double sx,sy,sz,red,green,blue,ireal;
	unsigned int si=0;

	// Read in spin coordinates and output to povray file
	for(int i=0; i<n_local_atoms;i++){
		spinfile >> sx >> sy >> sz;
        spinx[si]=sx;
        spiny[si]=sy;
        spinz[si]=sz;
		si++;

	}

		// close master file
	spinfile.close();

	// now read subsidiary files
	for(int file=0; file<n_files; file++){
		std::ifstream infile;
		infile.open(filenames[file].c_str());

		// read number of atoms in this file
		getline(infile,dummy);
		n_local_atoms=atoi(dummy.c_str());
		for(int i=0;i<n_local_atoms;i++){
			infile >> sx >> sy >> sz;
            spinx[si]=sx;
            spiny[si]=sy;
            spinz[si]=sz;
			si++;

		}
		// close subsidiary file
		infile.close();
	}

   // write .vtu file
   vtkfile << "<?xml version=\"1.0\"?>" << "\n";
   vtkfile << "<VTKFile type=\"UnstructuredGrid\">" << "\n";
   vtkfile << "<UnstructuredGrid>" << "\n";
   vtkfile << "<Piece NumberOfPoints=\""<<si<<"\"  NumberOfCells=\"1\">" << "\n";
   vtkfile << "<PointData Scalar=\"Spin\">" << "\n";
   vtkfile << "<DataArray type=\"Float32\" Name=\"Spin\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n";
   for(int i=0; i<si; ++i){
      vtkfile << spinx[i] << "\t" << spiny[i] << "\t" << spinz[i] << "\n";
   }
   vtkfile << "</DataArray>" << "\n";
   vtkfile << "</PointData>" << "\n";
   vtkfile << "<CellData>" << "\n";
   vtkfile << "</CellData>" << "\n";
   vtkfile << "<Points>" << "\n";
   vtkfile << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n";
   for(int i=0; i<si; ++i){
      vtkfile << cx[i] << "\t" << cy[i] << "\t" << cz[i] << "\n";
   }
   vtkfile << "</DataArray>" << "\n";
   vtkfile << "</Points>" << "\n";
   vtkfile << "<Cells>" << "\n";
   vtkfile << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << "\n";
   vtkfile << "1" << "\n";
   vtkfile << "</DataArray>" << "\n";
   vtkfile << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << "\n";
   vtkfile << "1" << "\n";
   vtkfile << "</DataArray>" << "\n";
   vtkfile << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << "\n";
   vtkfile << "1" << "\n";
   vtkfile << "</DataArray>" << "\n";
   vtkfile << "</Cells>" << "\n";
   vtkfile << "</Piece>" << "\n";
   vtkfile << "</UnstructuredGrid>" << "\n";
   vtkfile << "</VTKFile>" << "\n";

	// close vtk file
	vtkfile.close();

	spinfile_counter++;
      }
	  else ios=1;
   }



	// finished
	return 0;

}
