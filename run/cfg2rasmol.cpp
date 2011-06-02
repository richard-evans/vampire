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
	
	//std::vector <int> mat(0);
	//std::vector <int> cat(0); 
	//std::vector <double> cx(0);
	//std::vector <double> cy(0);
	//std::vector <double> cz(0);
	//std::vector <std::string> type(0);
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
	//std::vector <int> mat.resize(n_atoms);
	//std::vector <int> cat.resize(n_atoms); 
	//std::vector <double> cx.resize(n_atoms);
	//std::vector <double> cy.resize(n_atoms);
	//std::vector <double> cz.resize(n_atoms);
	//std::vector <std::string> type.resize(n_atoms);
	
	// Open output file
	std::ofstream outfile;
	outfile.open("crystal.xyz");
	
	// Output rasmol file header
	outfile << n_atoms << std::endl << std::endl;
	
	double cx,cy,cz;
	int mat, cat;
	
	std::cout << "here1" << std::endl;
	std::string atom_type;
	std::cout << "here2" << std::endl;
	
	std::cout << "Outputting data to rasmol file" << std::endl;
	
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
	//std::cout << n_atoms << "\t" << n_files << std::endl;
	/*
			// Set local output filename
		std::stringstream file_sstr;
		file_sstr << "atoms-coords";
		// Set CPUID on non-root process
		if(vmpi::my_rank!=0){
			file_sstr << "-" << std::setfill('0') << std::setw(5) << vmpi::my_rank;
		}
		file_sstr << ".cfg";
		std::string cfg_file = file_sstr.str();
		const char* cfg_filec = cfg_file.c_str();

		// Declare and open output file
		std::ofstream cfg_file_ofstr;
		cfg_file_ofstr.open (cfg_filec);
		
		// Output masterfile header on root process
		if(vmpi::my_rank==0){
			// Get system date
			char *asctime( const struct tm *time_ptr );

			cfg_file_ofstr << "#------------------------------------------------------"<< std::endl;
			cfg_file_ofstr << "# Atomistic coordinates configuration file for vampire"<< std::endl;
			cfg_file_ofstr << "#------------------------------------------------------"<< std::endl;
			cfg_file_ofstr << "# Date: "<< asctime << std::endl;
			cfg_file_ofstr << "#------------------------------------------------------"<< std::endl;
			cfg_file_ofstr << "Number of atoms: "<< vout::total_output_atoms << std::endl;
			cfg_file_ofstr << "#------------------------------------------------------" << std::endl;
			cfg_file_ofstr << "Number of spin files: " << vmpi::num_processors-1 << std::endl;
			for(int p=1;p<vmpi::num_processors;p++){
				std::stringstream cfg_sstr;
				cfg_sstr << "atoms-coords-" << std::setfill('0') << std::setw(5) << p << ".cfg";
				cfg_file_ofstr << cfg_sstr.str() << "\"" << std::endl;
			}
			cfg_file_ofstr << "#------------------------------------------------------"<< std::endl;
		}

		// Everyone now outputs their atom list
		cfg_file_ofstr << vout::local_output_atom_list.size() << std::endl;
	  	for(int i=0; i<vout::local_output_atom_list.size(); i++){
			const int atom = vout::local_output_atom_list[i];
			cfg_file_ofstr << atoms::type_array[atom] << "\t" << atoms::category_array[atom] << "\t" << 
			atoms::x_coord_array[atom] << "\t" << atoms::y_coord_array[atom] << "\t" << atoms::z_coord_array[atom] << "\t" << 
			mp::material[atoms::type_array[atom]].element << std::endl;
		}
	
		cfg_file_ofstr.close();*/
	
	// resize arrays
	
	// read in process 0 data
	
	// open and read in other process data
	
	// open and output crystal.xyz
	
}
