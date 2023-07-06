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
#include "exchange.hpp"
#include "material.hpp"
#include "unitcell.hpp"
#include "vio.hpp"

// unitcell module headers
#include "internal.hpp"

namespace unitcell{
namespace internal{

void read_unit_cell(unit_cell_t & unit_cell, std::string filename){

	std::cout << "Reading in unit cell data from disk..." << std::flush;
	zlog << zTs() << "Reading in unit cell data from disk..." << std::endl;

	// stringstream stream declaration
	std::stringstream inputfile;

   // fill input file stream with contents of file opened on master process
   inputfile.str( vin::get_string(filename.c_str(), "input", -1) );

   std::cout << "done!\nProcessing unit cell data..." << std::flush;
   zlog << zTs() << "Reading data completed. Processing unit cell data..." << std::endl;

	// keep record of current line
	unsigned int line_counter=0;
	unsigned int line_id=0;

   std::string exchange_type_string; // string defining exchange type

   // defaults for interaction list
   unsigned int interaction_range = 1; // assume +-1 unit cell as default

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

            std::cout << "\nProcessing data for " << unit_cell.atom.size() << " atoms..." << std::flush;
            zlog << zTs() << "\t" << "Processing data for " << unit_cell.atom.size() << " unit cell atoms..." << std::endl;


            // loop over all atoms and read into class
            for(unsigned int i = 0; i < unit_cell.atom.size(); i++){

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
                       << " of unit cell input file " << filename.c_str() << " is greater than the number of materials ( " << mp::num_materials << " ) specified in the material file. Exiting" << std::endl;
                  err::vexit();
               }
					unit_cell.atom[i].lc=lcat_id;
					unit_cell.atom[i].hc=hcat_id;
					//std::cout << i << "\t" << id << "\t" << cx << "\t" << cy << "\t" << cz << "\t" << mat_id << "\t" << lcat_id << "\t" << hcat_id << std::endl;
				}
				break;
			case 5:{

            // read (bilinear) exchange interactions
            unit_cell.bilinear.read_interactions(num_uc_atoms, inputfile, iss, filename, line_counter, interaction_range);
				break;

         }
         case 6:{

            // read biquadratic exchange interactions
            unit_cell.biquadratic.read_interactions(num_uc_atoms, inputfile, iss, filename, line_counter, interaction_range);
            break;

         }

			default:
				terminaltextcolor(RED);
				std::cerr << "Error! Unknown line type on line " << line_counter
					<< " of unit cell input file " << filename.c_str() << ". Exiting" << std::endl; err::vexit();
				terminaltextcolor(WHITE);
		}
		line_id++;
	} // end of while loop

   std::cout << "done!\nVerifying exchange interactions..." << std::flush;
   zlog << zTs() << "\t" << "Processing unit cell interactions completed" << std::endl;
   zlog << zTs() << "\t" << "Verifying unit cell exchange interactions..." << std::endl;

   // Verify exchange interactions are symmetric (required for MPI parallelization)
   unit_cell.bilinear.verify(filename);
   unit_cell.biquadratic.verify(filename);

   // If biquadratic interactins are included, then set flag to enable them
   if(unit_cell.biquadratic.interaction.size()>0){
      zlog << zTs() << "Enabling biquadratic interactions from unit cell file" << std::endl;
      exchange::biquadratic = true;
   }

   // set interaction range if larger than existing range
   if(interaction_range > unit_cell.interaction_range) unit_cell.interaction_range = interaction_range;

   std::cout << "done!" << std::endl;
   zlog << zTs() << "Verifying unit cell exchange interactions completed" << std::endl;
	zlog << zTs() << "\t" << "Number of atoms read-in: " << unit_cell.atom.size() << std::endl;
	zlog << zTs() << "\t" << "Number of bilinear interactions read-in: " << unit_cell.bilinear.interaction.size() << std::endl;
   zlog << zTs() << "\t" << "Number of biquadratic interactions read-in: " << unit_cell.biquadratic.interaction.size() << std::endl;
	zlog << zTs() << "\t" << "Exchange type: " << exchange_type_string << std::endl;
	zlog << zTs() << "\t" << "Calculated interaction range: " << unit_cell.interaction_range << " Unit Cells" << std::endl;

	return;
}

} // end of internal namespace
} // end of unitcell namespace
