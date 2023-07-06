//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard Evans 2022. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

// vampire headers
#include "create.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

// Internal create header
#include "internal.hpp"

namespace cs{

int create_crystal_structure(std::vector<cs::catom_t> & catom_array){
	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "cs::create_crystal_structure has been called" << std::endl;}

	int min_bounds[3];
	int max_bounds[3];

	#ifdef MPICF
	if(vmpi::mpi_mode==0){
		min_bounds[0] = int(vmpi::min_dimensions[0]/unit_cell.dimensions[0]);
		min_bounds[1] = int(vmpi::min_dimensions[1]/unit_cell.dimensions[1]);
		min_bounds[2] = int(vmpi::min_dimensions[2]/unit_cell.dimensions[2]);
		max_bounds[0] = vmath::iceil(vmpi::max_dimensions[0]/unit_cell.dimensions[0]);
		max_bounds[1] = vmath::iceil(vmpi::max_dimensions[1]/unit_cell.dimensions[1]);
		max_bounds[2] = vmath::iceil(vmpi::max_dimensions[2]/unit_cell.dimensions[2]);
	}
	else{
		min_bounds[0]=0;
		min_bounds[1]=0;
		min_bounds[2]=0;
		max_bounds[0]=cs::total_num_unit_cells[0];
		max_bounds[1]=cs::total_num_unit_cells[1];
		max_bounds[2]=cs::total_num_unit_cells[2];
	}
	#else
		min_bounds[0]=0;
		min_bounds[1]=0;
		min_bounds[2]=0;
		max_bounds[0]=cs::total_num_unit_cells[0];
		max_bounds[1]=cs::total_num_unit_cells[1];
		max_bounds[2]=cs::total_num_unit_cells[2];
	#endif

	cs::local_num_unit_cells[0]=max_bounds[0]-min_bounds[0];
	cs::local_num_unit_cells[1]=max_bounds[1]-min_bounds[1];
	cs::local_num_unit_cells[2]=max_bounds[2]-min_bounds[2];

	int num_atoms=cs::local_num_unit_cells[0]*cs::local_num_unit_cells[1]*cs::local_num_unit_cells[2]*unit_cell.atom.size();

	// set catom_array size
	catom_array.reserve(num_atoms);

	// Initialise atoms number
	int atom=0;

   // find maximum height lh_category
   unsigned int maxlh=0;
   for(unsigned int uca=0;uca<unit_cell.atom.size();uca++) if(unit_cell.atom[uca].hc > maxlh) maxlh = unit_cell.atom[uca].hc;
   maxlh+=1;

	// specify which atoms in the unit cell to generate based on whether they are defined in the material file
	std::vector<bool> inc_uc_atom(mp::max_materials, false);
	for( auto m : create::internal::mp) inc_uc_atom[m.unit_cell_category] = true;

	// Duplicate unit cell
	for(int z=min_bounds[2];z<max_bounds[2];z++){
		for(int y=min_bounds[1];y<max_bounds[1];y++){
			for(int x=min_bounds[0];x<max_bounds[0];x++){

				// need to change this to accept non-orthogonal lattices
				// Loop over atoms in unit cell
				for(unsigned int uca=0;uca<unit_cell.atom.size();uca++){
					double cx = (double(x)+unit_cell.atom[uca].x)*unit_cell.dimensions[0];
					double cy = (double(y)+unit_cell.atom[uca].y)*unit_cell.dimensions[1];
					double cz = (double(z)+unit_cell.atom[uca].z)*unit_cell.dimensions[2];
					#ifdef MPICF
						if(vmpi::mpi_mode==0){
							// only generate atoms within allowed dimensions
                     if(   (cx>=vmpi::min_dimensions[0] && cx<vmpi::max_dimensions[0]) &&
                           (cy>=vmpi::min_dimensions[1] && cy<vmpi::max_dimensions[1]) &&
                           (cz>=vmpi::min_dimensions[2] && cz<vmpi::max_dimensions[2])){
						#endif
							if(inc_uc_atom[unit_cell.atom[uca].mat] &&
								cx < cs::system_dimensions[0] &&
								cy < cs::system_dimensions[1] &&
								cz < cs::system_dimensions[2]
							){
							catom_array.push_back(cs::catom_t());
							catom_array[atom].x=cx;
							catom_array[atom].y=cy;
							catom_array[atom].z=cz;
							//std::cout << atom << "\t" << cx << "\t" << cy <<"\t" << cz << std::endl;
							catom_array[atom].material=unit_cell.atom[uca].mat;
							catom_array[atom].uc_id=uca;
							catom_array[atom].lh_category=unit_cell.atom[uca].hc+z*maxlh;
							catom_array[atom].uc_category=unit_cell.atom[uca].mat; // determine initial material (uc_category) for unit cell
							catom_array[atom].scx=x;
							catom_array[atom].scy=y;
							catom_array[atom].scz=z;
							atom++;
							}
						#ifdef MPICF
							}
						}
						else{
							if((cx<cs::system_dimensions[0]) && (cy<cs::system_dimensions[1]) && (cz<cs::system_dimensions[2])){
							catom_array.push_back(cs::catom_t());
							catom_array[atom].x=cx;
							catom_array[atom].y=cy;
							catom_array[atom].z=cz;
							catom_array[atom].material=unit_cell.atom[uca].mat;
							catom_array[atom].uc_id=uca;
							catom_array[atom].lh_category=unit_cell.atom[uca].hc+z*maxlh;
							catom_array[atom].uc_category=unit_cell.atom[uca].mat; // determine initial material (uc_category) for unit cell
							catom_array[atom].scx=x;
							catom_array[atom].scy=y;
							catom_array[atom].scz=z;
							catom_array[atom].include=false; // assume no atoms until classification complete
							atom++;
							}
						}
						#endif
					}
				}
			}
		}

	// Check to see if actual and expected number of atoms agree, if not trim the excess
	if(atom!=num_atoms){
		std::vector<cs::catom_t> tmp_catom_array(num_atoms);
		tmp_catom_array=catom_array;
		catom_array.resize(atom);
		for(int a=0;a<atom;a++){
			catom_array[a]=tmp_catom_array[a];
		}
		tmp_catom_array.resize(0);
	}

	// Check to see if any atoms have been generated
	if(atom==0){
		terminaltextcolor(RED);
		std::cout << "Error - no atoms have been generated, increase system dimensions!" << std::endl;
		terminaltextcolor(WHITE);
		zlog << zTs() << "Error: No atoms have been generated. Increase system dimensions." << std::endl;
		err::vexit();
	}

	// Now unselect all atoms by default for particle shape cutting
	for(unsigned int atom=0;atom<catom_array.size();atom++){
		catom_array[atom].include=false;
	}

	return EXIT_SUCCESS;
}







} // end of namespace cs
