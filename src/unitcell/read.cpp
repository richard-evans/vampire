//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <sstream>

// Vampire headers
#include "errors.hpp"
#include "material.hpp"
#include "unitcell.hpp"
#include "vio.hpp"

// unitcell module headers
#include "internal.hpp"

namespace unitcell{
namespace internal{

void read_unit_cell(unit_cell_t & unit_cell, std::string filename){

	std::cout << "Reading in unit cell data..." << std::flush;
	zlog << zTs() << "Reading in unit cell data..." << std::endl;

	// stringstream stream declaration
	std::stringstream inputfile;

   // fill input file stream with contents of file opened on master process
   inputfile.str( vin::get_string(filename.c_str(), "input", -1) );

	// keep record of current line
	unsigned int line_counter=0;
	unsigned int line_id=0;
	// Loop over all lines
	while (! inputfile.eof() ){
		line_counter++;
		// read in whole line
		std::string line;
		getline(inputfile,line);
		//std::cout << line.c_str() << std::endl;

		// ignore blank lines
		std::string empty="";
		if(line==empty) continue;

		// set character triggers
		const char* hash="#";	// Comment identifier

		bool has_hash=false;
		// Determine if line is a comment line
		for(unsigned int i=0;i<line.length();i++){
			char c=line.at(i);

			if(c== *hash){
					has_hash=true;
					break;
			}
		}
		// if hash character found then read next line
		if(has_hash==true) continue;

		// convert line to string stream
		std::istringstream iss(line,std::istringstream::in);

		// defaults for interaction list
		int exc_type=-1; // assume isotropic
		int num_interactions=0; // assume no interactions
		int interaction_range=1; // assume +-1 unit cell as default

		// non-comment line found - check for line number
		switch(line_id){
			case 0:
				iss >> unit_cell.dimensions[0] >> unit_cell.dimensions[1] >> unit_cell.dimensions[2];
				break;
			case 1:
				iss >> unit_cell.shape[0][0] >> unit_cell.shape[0][1] >> unit_cell.shape[0][2];
				break;
			case 2:
				iss >> unit_cell.shape[1][0] >> unit_cell.shape[1][1] >> unit_cell.shape[1][2];
				break;
			case 3:
				iss >> unit_cell.shape[2][0] >> unit_cell.shape[2][1] >> unit_cell.shape[2][2];
				break;
			case 4:
				int num_uc_atoms;
				iss >> num_uc_atoms;
				//std::cout << "Reading in " << num_uc_atoms << " atoms" << std::endl;
				// resize unit_cell.atom array if within allowable bounds
				if( (num_uc_atoms >0) && (num_uc_atoms <= 1000000)) unit_cell.atom.resize(num_uc_atoms);
				else {
					terminaltextcolor(RED);
					std::cerr << "Error! Requested number of atoms " << num_uc_atoms << " on line " << line_counter
					<< " of unit cell input file " << filename.c_str() << " is outside of valid range 1-1,000,000. Exiting" << std::endl; err::vexit();
					terminaltextcolor(WHITE);
				}
				// loop over all atoms and read into class
				for (unsigned int i=0; i<unit_cell.atom.size(); i++){
					line_counter++;
					// declare safe temporaries for atom input
					int id=i;
					double cx=2.0, cy=2.0,cz=2.0; // coordinates - default will give an error
					int mat_id=0, lcat_id=0, hcat_id=0; // sensible defaults if omitted
					// get line
					std::string atom_line;
					getline(inputfile,atom_line);
					std::istringstream atom_iss(atom_line,std::istringstream::in);
					atom_iss >> id >> cx >> cy >> cz >> mat_id >> lcat_id >> hcat_id;
					//std::cout << id << "\t" << cx << "\t" << cy << "\t" << cz<< "\t"  << mat_id << "\t" << lcat_id << "\t" << hcat_id << std::endl;
					//inputfile >> id >> cx >> cy >> cz >> mat_id >> lcat_id >> hcat_id;
					// now check for mostly sane input
					if(cx>=0.0 && cx <=1.0) unit_cell.atom[i].x=cx;
					else{
						terminaltextcolor(RED);
						std::cerr << "Error! atom x-coordinate for atom " << id << " on line " << line_counter
									 << " of unit cell input file " << filename.c_str() << " is outside of valid range 0.0-1.0. Exiting" << std::endl;
						terminaltextcolor(WHITE);
						zlog << zTs() << "Error! atom x-coordinate for atom " << id << " on line " << line_counter
									 << " of unit cell input file " << filename.c_str() << " is outside of valid range 0.0-1.0. Exiting" << std::endl;
						err::vexit();
					}
					if(cy>=0.0 && cy <=1.0) unit_cell.atom[i].y=cy;
					else{
						terminaltextcolor(RED);
						std::cerr << "Error! atom y-coordinate for atom " << id << " on line " << line_counter
									 << " of unit cell input file " << filename.c_str() << " is outside of valid range 0.0-1.0. Exiting" << std::endl;
						terminaltextcolor(WHITE);
						zlog << zTs() << "Error! atom y-coordinate for atom " << id << " on line " << line_counter
									     << " of unit cell input file " << filename.c_str() << " is outside of valid range 0.0-1.0. Exiting" << std::endl;
						err::vexit();
					}
					if(cz>=0.0 && cz <=1.0) unit_cell.atom[i].z=cz;
					else{
						terminaltextcolor(RED);
						std::cerr << "Error! atom z-coordinate for atom " << id << " on line " << line_counter
						<< " of unit cell input file " << filename.c_str() << " is outside of valid range 0.0-1.0. Exiting" << std::endl;
						terminaltextcolor(WHITE);
						zlog << zTs() << "Error! atom z-coordinate for atom " << id << " on line " << line_counter
										  << " of unit cell input file " << filename.c_str() << " is outside of valid range 0.0-1.0. Exiting" << std::endl;
						err::vexit();
					}
					if(mat_id >=0 && mat_id<mp::num_materials) unit_cell.atom[i].mat=mat_id;
					else{
						terminaltextcolor(RED);
						std::cerr << "Error! Requested material id " << mat_id << " for atom number " << id <<  " on line " << line_counter
									 << " of unit cell input file " << filename.c_str() << " is greater than the number of materials ( " << mp::num_materials << " ) specified in the material file. Exiting" << std::endl;
						terminaltextcolor(WHITE);
						zlog << zTs() << "Error! Requested material id " << mat_id << " for atom number " << id <<  " on line " << line_counter
                            << " of unit cell input file " << filename.c_str() << " is greater than the number of materials ( " << mp::num_materials << " ) specified in the material file. Exiting" << std::endl; err::vexit();}
					unit_cell.atom[i].lc=lcat_id;
					unit_cell.atom[i].hc=hcat_id;
					//std::cout << i << "\t" << id << "\t" << cx << "\t" << cy << "\t" << cz << "\t" << mat_id << "\t" << lcat_id << "\t" << hcat_id << std::endl;
				}
				break;
			case 5:
				iss >> num_interactions >> exc_type;
				//std::cout << num_interactions << "\t" << exc_type << std::endl;
				if(num_interactions>=0) unit_cell.interaction.resize(num_interactions);
				else {
					terminaltextcolor(RED);
					std::cerr << "Error! Requested number of interactions " << num_interactions << " on line " << line_counter
					<< " of unit cell input file " << filename.c_str() << " is less than 0. Exiting" << std::endl; err::vexit();
				    terminaltextcolor(WHITE);
				}
				// if exchange type omitted, then assume isotropic values from material file
				//if(exc_type==-1) unit_cell.exchange_type=0;
				// loop over all interactions and read into class
				for (int i=0; i<num_interactions; i++){
					//std::cout << "setting up interaction "<< i+1<< " of " << num_interactions << " interactions" << std::endl;
					// declare safe temporaries for interaction input
					int id=i;
					int iatom=-1,jatom=-1; // atom pairs
					int dx=0, dy=0,dz=0; // relative unit cell coordinates
					// get line
					std::string int_line;
					getline(inputfile,int_line);
					//std::cout << int_line.c_str() << std::endl;
					std::istringstream int_iss(int_line,std::istringstream::in);
					int_iss >> id >> iatom >> jatom >> dx >> dy >> dz;
					//inputfile >> id >> iatom >> jatom >> dx >> dy >> dz;
					line_counter++;
					// check for sane input
					if(iatom>=0 && iatom < int(unit_cell.atom.size())) unit_cell.interaction[i].i=iatom;
					else if(iatom>=0 && iatom >= int(unit_cell.atom.size())){
						terminaltextcolor(RED);
						std::cerr << std::endl << "Error! iatom number "<< iatom <<" for interaction id " << id << " on line " << line_counter
							  << " of unit cell input file " << filename.c_str() << " is outside of valid range 0-"
							  << unit_cell.atom.size()-1 << ". Exiting" << std::endl;
						terminaltextcolor(WHITE);
						zlog << zTs() << "Error! iatom number "<< iatom <<" for interaction id " << id << " on line " << line_counter
						     << " of unit cell input file " << filename.c_str() << " is outside of valid range 0-"<< unit_cell.atom.size()-1
						     << ". Exiting" << std::endl;
						err::vexit();
					}
					else{
					  terminaltextcolor(RED);
					  std::cerr << std::endl << "Error! No valid interaction for interaction id " << id << " on line " << line_counter
						    << " of unit cell input file " << filename.c_str() << ". Possibly too many interactions defined. Exiting" << std::endl;
					  terminaltextcolor(WHITE);
					  zlog << zTs() << "Error! No valid interaction for interaction id " << id << " on line " << line_counter
					       << " of unit cell input file " << filename.c_str() << ". Possibly too many interactions defined. Exiting" << std::endl;
					  err::vexit();
					}
					if(iatom>=0 && jatom < int(unit_cell.atom.size())) unit_cell.interaction[i].j=jatom;
					else{
						terminaltextcolor(RED);
						std::cerr << std::endl << "Error! jatom number "<< jatom <<" for interaction id " << id << " on line " << line_counter
							  << " of unit cell input file " << filename.c_str() << " is outside of valid range 0-"
							  << unit_cell.atom.size()-1 << ". Exiting" << std::endl;
						terminaltextcolor(WHITE);
						zlog << zTs() << "Error! jatom number "<< jatom <<" for interaction id " << id << " on line " << line_counter
							  << " of unit cell input file " << filename.c_str() << " is outside of valid range 0-"
							  << unit_cell.atom.size()-1 << ". Exiting" << std::endl;
						err::vexit();
						}
					unit_cell.interaction[i].dx=dx;
					unit_cell.interaction[i].dy=dy;
					unit_cell.interaction[i].dz=dz;
					// check for long range interactions
					if(abs(dx)>interaction_range) interaction_range=abs(dx);
					if(abs(dy)>interaction_range) interaction_range=abs(dy);
					if(abs(dz)>interaction_range) interaction_range=abs(dz);

					int iatom_mat = unit_cell.atom[iatom].mat;
					int jatom_mat = unit_cell.atom[jatom].mat;
					switch(exc_type){
						//case -1: // assume isotropic
						//	unit_cell.interaction[i].Jij[0][0]=mp::material[iatom_mat].Jij_matrix[jatom_mat][0]; // only works if read after mat file
						//	break;
						case 0:
							int_iss >> unit_cell.interaction[i].Jij[0][0];
							//std::cout << i << "\t" << unit_cell.interaction[i].Jij[0][0] << std::endl;
							break;
						case 1:
							int_iss >> unit_cell.interaction[i].Jij[0][0] >> unit_cell.interaction[i].Jij[1][1] >> unit_cell.interaction[i].Jij[2][2];
							break;
						case 2:
							int_iss >> unit_cell.interaction[i].Jij[0][0] >> unit_cell.interaction[i].Jij[0][1] >> unit_cell.interaction[i].Jij[0][2];
							int_iss >> unit_cell.interaction[i].Jij[1][0] >> unit_cell.interaction[i].Jij[1][1] >> unit_cell.interaction[i].Jij[1][2];
							int_iss >> unit_cell.interaction[i].Jij[2][0] >> unit_cell.interaction[i].Jij[2][1] >> unit_cell.interaction[i].Jij[2][2];
							break;
						default:
							terminaltextcolor(RED);
							std::cerr << "Error! Requested exchange type " << exc_type << " on line " << line_counter
					<< " of unit cell input file " << filename.c_str() << " is outside of valid range 0-2. Exiting" << std::endl; err::vexit();
							terminaltextcolor(WHITE);
					}
					// increment number of interactions for atom i
					unit_cell.atom[iatom].ni++;
				}
				// set interaction range
				unit_cell.interaction_range=interaction_range;
				// set exchange type
				unit_cell.exchange_type=exc_type;
				break;
			default:
				terminaltextcolor(RED);
				std::cerr << "Error! Unknown line type on line " << line_counter
					<< " of unit cell input file " << filename.c_str() << ". Exiting" << std::endl; err::vexit();
				terminaltextcolor(WHITE);
		}
		line_id++;
	} // end of while loop

   // Verify exchange interactions are symmetric (required for MPI parallelization)
   uc::internal::verify_exchange_interactions(unit_cell, filename);

   std::cout << "Done!" << std::endl;
   zlog << "Done!" << std::endl;
	zlog << zTs() << "\t" << "Number of atoms read-in: " << unit_cell.atom.size() << std::endl;
	zlog << zTs() << "\t" << "Number of interactions read-in: " << unit_cell.interaction.size() << std::endl;
	zlog << zTs() << "\t" << "Exchange type: " <<  unit_cell.exchange_type << std::endl;
	zlog << zTs() << "\t" << "Calculated interaction range: " << unit_cell.interaction_range << " Unit Cells" << std::endl;

	return;
}

} // end of internal namespace
} // end of unitcell namespace
