/// Program to convert vampire cfg files to rasmol format
///
/// ./cfg2rasmol 

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


int main(){
	
	std::vector <std::string> filenames(0);
	
	// open coordinate file
	std::ifstream coord_file;
	coord_file.open("atoms-coords.cfg");
	
	// check for open file
	if(!coord_file.is_open()){
		std::cerr << "Error! Coordinate file atoms-coords.cfg cannot be opened. Exiting" << std::endl;
      exit(1);
	}
	
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
	
	// Open output file
	std::ofstream outfile;
	outfile.open("crystal.xyz");
	
	// Output rasmol file header
	outfile << n_atoms << std::endl << std::endl;
	
	double cx,cy,cz;
	int mat, cat;
	
	std::string atom_type;
	
	//std::cout << "Outputting data to rasmol file" << std::endl;
	
	// finish reading master file coordinates
	for(int i=0;i<n_local_atoms;i++){
		coord_file >> mat >> cat >> cx >> cy >> cz >> atom_type;
		outfile << atom_type << "\t" << cx << "\t" << cy << "\t" << cz << std::endl;
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
			infile >> mat >> cat >> cx >> cy >> cz >> atom_type;
			outfile << atom_type << "\t" << cx << "\t" << cy << "\t" << cz << std::endl;
		}
		// close subsidiary file
		infile.close();
	}
	
	// close rasmol file
	outfile.close();
	
	// finished
	return 0;

}
