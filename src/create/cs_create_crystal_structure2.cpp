///
/// @file 
/// @brief Crystal Structure generation routine
///
/// Replicates unit cell to fill system dimensions and adds various atomic identifiers
/// such as material, category
///
/// @author Richard Evans, richard.evans@york.ac.uk
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans 2009-2011. All Rights Reserved.
///
/// @internal
/// Created  08/06/2009
/// Revision  3.0
/// Copyright  Copyright (c) 2011, Richard Evans
///
///=====================================================================================
///

// Vampire Header files
#include "create.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

// Standard Libraries
#include <string>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>

namespace cs{
	
void unit_cell_set(cs::unit_cell_t &);
	
int create_crystal_structure(std::vector<cs::catom_t> & catom_array){
	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "cs::create_crystal_structure has been called" << std::endl;}	

	// populate unit cell coordinates
	unit_cell_set(cs::unit_cell);

	// Calculate number of global and local unit cells required (rounding up)
	cs::total_num_unit_cells[0]=int(vmath::iceil(cs::system_dimensions[0]/cs::unit_cell_size[0]));
	cs::total_num_unit_cells[1]=int(vmath::iceil(cs::system_dimensions[1]/cs::unit_cell_size[1]));
	cs::total_num_unit_cells[2]=int(vmath::iceil(cs::system_dimensions[2]/cs::unit_cell_size[2]));
	
	int min_bounds[3];
	int max_bounds[3];
	
	#ifdef MPICF
	if(vmpi::mpi_mode==0){
		min_bounds[0] = int(vmpi::min_dimensions[0]/cs::unit_cell_size[0]);
		min_bounds[1] = int(vmpi::min_dimensions[1]/cs::unit_cell_size[1]);
		min_bounds[2] = int(vmpi::min_dimensions[2]/cs::unit_cell_size[2]);
		max_bounds[0] = vmath::iceil(vmpi::max_dimensions[0]/cs::unit_cell_size[0]);
		max_bounds[1] = vmath::iceil(vmpi::max_dimensions[1]/cs::unit_cell_size[1]);
		max_bounds[2] = vmath::iceil(vmpi::max_dimensions[2]/cs::unit_cell_size[2]);
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

	unsigned int maxlh=2;

	// Duplicate unit cell
	for(int z=min_bounds[2];z<max_bounds[2];z++){
		for(int y=min_bounds[1];y<max_bounds[1];y++){
			for(int x=min_bounds[0];x<max_bounds[0];x++){
				// Loop over atoms in unit cell
				for(int uca=0;uca<unit_cell.atom.size();uca++){
					double cx = (double(x)+unit_cell.atom[uca].x)*cs::unit_cell_size[0];
					double cy = (double(y)+unit_cell.atom[uca].y)*cs::unit_cell_size[1];
					double cz = (double(z)+unit_cell.atom[uca].z)*cs::unit_cell_size[2];

					#ifdef MPICF
						if(vmpi::mpi_mode==0){
							// only generate atoms within allowed dimensions
							if(	(cx>=vmpi::min_dimensions[0] && cx<vmpi::max_dimensions[0]) &&
									(cy>=vmpi::min_dimensions[1] && cy<vmpi::max_dimensions[1]) &&
									(cz>=vmpi::min_dimensions[2] && cz<vmpi::max_dimensions[2])){
						#endif
							if((cx<=cs::system_dimensions[0]) && (cy<=cs::system_dimensions[1]) && (cz<=cs::system_dimensions[2])){
							catom_array.push_back(cs::catom_t());
							catom_array[atom].x=cx;
							catom_array[atom].y=cy;
							catom_array[atom].z=cz;
							//std::cout << atom << "\t" << cx << "\t" << cy <<"\t" << cz << std::endl;
							catom_array[atom].material=0;
							catom_array[atom].lh_category=unit_cell.atom[uca].hc+z*maxlh;
							catom_array[atom].uc_category=unit_cell.atom[uca].lc;
							catom_array[atom].scx=x;
							catom_array[atom].scy=y;
							catom_array[atom].scz=z;
							atom++;
							}
						#ifdef MPICF
							}
						}
						else{
							if((cx<=cs::system_dimensions[0]) && (cy<=cs::system_dimensions[1]) && (cz<=cs::system_dimensions[2])){
							catom_array.push_back(cs::catom_t());
							catom_array[atom].x=cx;
							catom_array[atom].y=cy;
							catom_array[atom].z=cz;
							catom_array[atom].material=0;
							catom_array[atom].material=0;
							catom_array[atom].lh_category=unit_cell.atom[uca].hc+z*maxlh;
							catom_array[atom].uc_category=unit_cell.atom[uca].lc;
							catom_array[atom].scx=x;
							catom_array[atom].scy=y;
							catom_array[atom].scz=z;
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
	
	// determine z-bounds for materials
	std::vector<double> mat_min(mp::num_materials);
	std::vector<double> mat_max(mp::num_materials);

	for(int mat=0;mat<mp::num_materials;mat++){
		mat_min[mat]=mp::material[mat].min*cs::system_dimensions[2];
		mat_max[mat]=mp::material[mat].max*cs::system_dimensions[2];
		if(mat_max[mat]<0.0000001) mat_max[mat]=-0.1;
	}
	
	// Assign materials to generated atoms
	for(unsigned int atom=0;atom<catom_array.size();atom++){
		for(int mat=0;mat<mp::num_materials;mat++){
			const double cz=catom_array[atom].z;
			if((cz>=mat_min[mat]) && (cz<mat_max[mat])) catom_array[atom].material=mat;
		}
	}

	// Check to see if any atoms have been generated
	if(atom==0){
		std::cout << "Error - no atoms have been generated, increase system dimensions!" << std::endl;
		err::vexit();
	}

	return EXIT_SUCCESS;
}

void unit_cell_set(unit_cell_t & unit_cell){

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "cs::unit_cell_set has been called" << std::endl;}	

		// Simple Cubic
		if(cs::crystal_structure=="sc"){
			unit_cell.size=1;
			unit_cell.lcsize=1;
			unit_cell.hcsize=1;
			unit_cell.interaction_range=1;
			unit_cell.atom.resize(1);
			//-----------------------------
			unit_cell.atom[0].x=0.0;
			unit_cell.atom[0].y=0.0;
			unit_cell.atom[0].z=0.0;
			unit_cell.atom[0].lc=0;
			unit_cell.atom[0].hc=0;
			//-----------------------------
			unit_cell.interaction.resize(6);
			//-----------------------------
			unit_cell.interaction[0].i=0;
			unit_cell.interaction[0].j=0;
			unit_cell.interaction[0].dx=1;
			unit_cell.interaction[0].dy=0;
			unit_cell.interaction[0].dz=0;
			//-----------------------------
			unit_cell.interaction[1].i=0;
			unit_cell.interaction[1].j=0;
			unit_cell.interaction[1].dx=-1;
			unit_cell.interaction[1].dy=0;
			unit_cell.interaction[1].dz=0;
			//-----------------------------
			unit_cell.interaction[2].i=0;
			unit_cell.interaction[2].j=0;
			unit_cell.interaction[2].dx=0;
			unit_cell.interaction[2].dy=1;
			unit_cell.interaction[2].dz=0;
			//-----------------------------
			unit_cell.interaction[3].i=0;
			unit_cell.interaction[3].j=0;
			unit_cell.interaction[3].dx=0;
			unit_cell.interaction[3].dy=-1;
			unit_cell.interaction[3].dz=0;
			//-----------------------------
			unit_cell.interaction[4].i=0;
			unit_cell.interaction[4].j=0;
			unit_cell.interaction[4].dx=0;
			unit_cell.interaction[4].dy=0;
			unit_cell.interaction[4].dz=1;
			//-----------------------------
			unit_cell.interaction[5].i=0;
			unit_cell.interaction[5].j=0;
			unit_cell.interaction[5].dx=0;
			unit_cell.interaction[5].dy=0;
			unit_cell.interaction[5].dz=-1;
			//-----------------------------
		}
		else if(cs::crystal_structure=="bcc"){
			unit_cell.size=2;
			unit_cell.lcsize=2;
			unit_cell.hcsize=2;
			unit_cell.interaction_range=1;
			unit_cell.atom.resize(2);
			//-----------------------------
			unit_cell.atom[0].x=0.0;
			unit_cell.atom[0].y=0.0;
			unit_cell.atom[0].z=0.0;
			unit_cell.atom[0].lc=0;
			unit_cell.atom[0].hc=0;
			//-----------------------------
			unit_cell.atom[1].x=0.5;
			unit_cell.atom[1].y=0.5;
			unit_cell.atom[1].z=0.5;
			unit_cell.atom[1].lc=1;
			unit_cell.atom[1].hc=1;
			//-----------------------------
		}
		else if(cs::crystal_structure=="fct"){
			unit_cell.size=2;
			unit_cell.lcsize=2;
			unit_cell.hcsize=1;
			unit_cell.interaction_range=1;
			unit_cell.atom.resize(2);
			//-----------------------------
			unit_cell.atom[0].x=0.0;
			unit_cell.atom[0].y=0.0;
			unit_cell.atom[0].z=0.0;
			unit_cell.atom[0].lc=0;
			unit_cell.atom[0].hc=0;
			//-----------------------------
			unit_cell.atom[1].x=0.5;
			unit_cell.atom[1].y=0.5;
			unit_cell.atom[1].z=0.0;
			unit_cell.atom[1].lc=1;
			unit_cell.atom[1].hc=0;
			//-----------------------------
		}
		else if(cs::crystal_structure=="fcc"){
			unit_cell.size=4;
			unit_cell.lcsize=4;
			unit_cell.hcsize=2;
			unit_cell.interaction_range=1;
			unit_cell.atom.resize(4);
			//-----------------------------
			unit_cell.atom[0].x=0.0;
			unit_cell.atom[0].y=0.0;
			unit_cell.atom[0].z=0.0;
			unit_cell.atom[0].lc=0;
			unit_cell.atom[0].hc=0;
			//-----------------------------
			unit_cell.atom[1].x=0.5;
			unit_cell.atom[1].y=0.5;
			unit_cell.atom[1].z=0.0;
			unit_cell.atom[1].lc=1;
			unit_cell.atom[1].hc=0;
			//-----------------------------
			unit_cell.atom[2].x=0.5;
			unit_cell.atom[2].y=0.0;
			unit_cell.atom[2].z=0.5;
			unit_cell.atom[2].lc=2;
			unit_cell.atom[2].hc=1;
			//-----------------------------
			unit_cell.atom[3].x=0.0;
			unit_cell.atom[3].y=0.5;
			unit_cell.atom[3].z=0.5;
			unit_cell.atom[3].lc=3;
			unit_cell.atom[3].hc=1;
			//-----------------------------
			unit_cell.interaction.resize(48);
			//-----------------------------
			unit_cell.interaction[0].i=0;
			unit_cell.interaction[0].j=1;
			unit_cell.interaction[0].dx=-1;
			unit_cell.interaction[0].dy=-1;
			unit_cell.interaction[0].dz=0;
			//----------------------------------------------
			unit_cell.interaction[1].i=0;
			unit_cell.interaction[1].j=2;
			unit_cell.interaction[1].dx=-1;
			unit_cell.interaction[1].dy=0;
			unit_cell.interaction[1].dz=-1;
			//----------------------------------------------
			unit_cell.interaction[2].i=0;
			unit_cell.interaction[2].j=1;
			unit_cell.interaction[2].dx=-1;
			unit_cell.interaction[2].dy=0;
			unit_cell.interaction[2].dz=0;
			//----------------------------------------------
			unit_cell.interaction[3].i=0;
			unit_cell.interaction[3].j=2;
			unit_cell.interaction[3].dx=-1;
			unit_cell.interaction[3].dy=0;
			unit_cell.interaction[3].dz=0;
			//----------------------------------------------
			unit_cell.interaction[4].i=0;
			unit_cell.interaction[4].j=3;
			unit_cell.interaction[4].dx=0;
			unit_cell.interaction[4].dy=-1;
			unit_cell.interaction[4].dz=-1;
			//----------------------------------------------
			unit_cell.interaction[5].i=0;
			unit_cell.interaction[5].j=1;
			unit_cell.interaction[5].dx=0;
			unit_cell.interaction[5].dy=-1;
			unit_cell.interaction[5].dz=0;
			//----------------------------------------------
			unit_cell.interaction[6].i=0;
			unit_cell.interaction[6].j=3;
			unit_cell.interaction[6].dx=0;
			unit_cell.interaction[6].dy=-1;
			unit_cell.interaction[6].dz=0;
			//----------------------------------------------
			unit_cell.interaction[7].i=0;
			unit_cell.interaction[7].j=2;
			unit_cell.interaction[7].dx=0;
			unit_cell.interaction[7].dy=0;
			unit_cell.interaction[7].dz=-1;
			//----------------------------------------------
			unit_cell.interaction[8].i=0;
			unit_cell.interaction[8].j=3;
			unit_cell.interaction[8].dx=0;
			unit_cell.interaction[8].dy=0;
			unit_cell.interaction[8].dz=-1;
			//----------------------------------------------
			unit_cell.interaction[9].i=0;
			unit_cell.interaction[9].j=1;
			unit_cell.interaction[9].dx=0;
			unit_cell.interaction[9].dy=0;
			unit_cell.interaction[9].dz=0;
			//----------------------------------------------
			unit_cell.interaction[10].i=0;
			unit_cell.interaction[10].j=2;
			unit_cell.interaction[10].dx=0;
			unit_cell.interaction[10].dy=0;
			unit_cell.interaction[10].dz=0;
			//----------------------------------------------
			unit_cell.interaction[11].i=0;
			unit_cell.interaction[11].j=3;
			unit_cell.interaction[11].dx=0;
			unit_cell.interaction[11].dy=0;
			unit_cell.interaction[11].dz=0;
			//----------------------------------------------
			unit_cell.interaction[12].i=1;
			unit_cell.interaction[12].j=2;
			unit_cell.interaction[12].dx=0;
			unit_cell.interaction[12].dy=0;
			unit_cell.interaction[12].dz=-1;
			//----------------------------------------------
			unit_cell.interaction[13].i=1;
			unit_cell.interaction[13].j=3;
			unit_cell.interaction[13].dx=0;
			unit_cell.interaction[13].dy=0;
			unit_cell.interaction[13].dz=-1;
			//----------------------------------------------
			unit_cell.interaction[14].i=1;
			unit_cell.interaction[14].j=0;
			unit_cell.interaction[14].dx=0;
			unit_cell.interaction[14].dy=0;
			unit_cell.interaction[14].dz=0;
			//----------------------------------------------
			unit_cell.interaction[15].i=1;
			unit_cell.interaction[15].j=2;
			unit_cell.interaction[15].dx=0;
			unit_cell.interaction[15].dy=0;
			unit_cell.interaction[15].dz=0;
			//----------------------------------------------
			unit_cell.interaction[16].i=1;
			unit_cell.interaction[16].j=3;
			unit_cell.interaction[16].dx=0;
			unit_cell.interaction[16].dy=0;
			unit_cell.interaction[16].dz=0;
			//----------------------------------------------
			unit_cell.interaction[17].i=1;
			unit_cell.interaction[17].j=2;
			unit_cell.interaction[17].dx=0;
			unit_cell.interaction[17].dy=1;
			unit_cell.interaction[17].dz=-1;
			//----------------------------------------------
			unit_cell.interaction[18].i=1;
			unit_cell.interaction[18].j=0;
			unit_cell.interaction[18].dx=0;
			unit_cell.interaction[18].dy=1;
			unit_cell.interaction[18].dz=0;
			//----------------------------------------------
			unit_cell.interaction[19].i=1;
			unit_cell.interaction[19].j=2;
			unit_cell.interaction[19].dx=0;
			unit_cell.interaction[19].dy=1;
			unit_cell.interaction[19].dz=0;
			//----------------------------------------------
			unit_cell.interaction[20].i=1;
			unit_cell.interaction[20].j=3;
			unit_cell.interaction[20].dx=1;
			unit_cell.interaction[20].dy=0;
			unit_cell.interaction[20].dz=-1;
			//----------------------------------------------
			unit_cell.interaction[21].i=1;
			unit_cell.interaction[21].j=0;
			unit_cell.interaction[21].dx=1;
			unit_cell.interaction[21].dy=0;
			unit_cell.interaction[21].dz=0;
			//----------------------------------------------
			unit_cell.interaction[22].i=1;
			unit_cell.interaction[22].j=3;
			unit_cell.interaction[22].dx=1;
			unit_cell.interaction[22].dy=0;
			unit_cell.interaction[22].dz=0;
			//----------------------------------------------
			unit_cell.interaction[23].i=1;
			unit_cell.interaction[23].j=0;
			unit_cell.interaction[23].dx=1;
			unit_cell.interaction[23].dy=1;
			unit_cell.interaction[23].dz=0;
			//----------------------------------------------
			unit_cell.interaction[24].i=2;
			unit_cell.interaction[24].j=1;
			unit_cell.interaction[24].dx=0;
			unit_cell.interaction[24].dy=-1;
			unit_cell.interaction[24].dz=0;
			//----------------------------------------------
			unit_cell.interaction[25].i=2;
			unit_cell.interaction[25].j=3;
			unit_cell.interaction[25].dx=0;
			unit_cell.interaction[25].dy=-1;
			unit_cell.interaction[25].dz=0;
			//----------------------------------------------
			unit_cell.interaction[26].i=2;
			unit_cell.interaction[26].j=1;
			unit_cell.interaction[26].dx=0;
			unit_cell.interaction[26].dy=-1;
			unit_cell.interaction[26].dz=1;
			//----------------------------------------------
			unit_cell.interaction[27].i=2;
			unit_cell.interaction[27].j=0;
			unit_cell.interaction[27].dx=0;
			unit_cell.interaction[27].dy=0;
			unit_cell.interaction[27].dz=0;
			//----------------------------------------------
			unit_cell.interaction[28].i=2;
			unit_cell.interaction[28].j=1;
			unit_cell.interaction[28].dx=0;
			unit_cell.interaction[28].dy=0;
			unit_cell.interaction[28].dz=0;
			//----------------------------------------------
			unit_cell.interaction[29].i=2;
			unit_cell.interaction[29].j=3;
			unit_cell.interaction[29].dx=0;
			unit_cell.interaction[29].dy=0;
			unit_cell.interaction[29].dz=0;
			//----------------------------------------------
			unit_cell.interaction[30].i=2;
			unit_cell.interaction[30].j=0;
			unit_cell.interaction[30].dx=0;
			unit_cell.interaction[30].dy=0;
			unit_cell.interaction[30].dz=1;
			//----------------------------------------------
			unit_cell.interaction[31].i=2;
			unit_cell.interaction[31].j=1;
			unit_cell.interaction[31].dx=0;
			unit_cell.interaction[31].dy=0;
			unit_cell.interaction[31].dz=1;
			//----------------------------------------------
			unit_cell.interaction[32].i=2;
			unit_cell.interaction[32].j=3;
			unit_cell.interaction[32].dx=1;
			unit_cell.interaction[32].dy=-1;
			unit_cell.interaction[32].dz=0;
			//----------------------------------------------
			unit_cell.interaction[33].i=2;
			unit_cell.interaction[33].j=0;
			unit_cell.interaction[33].dx=1;
			unit_cell.interaction[33].dy=0;
			unit_cell.interaction[33].dz=0;
			//----------------------------------------------
			unit_cell.interaction[34].i=2;
			unit_cell.interaction[34].j=3;
			unit_cell.interaction[34].dx=1;
			unit_cell.interaction[34].dy=0;
			unit_cell.interaction[34].dz=0;
			//----------------------------------------------
			unit_cell.interaction[35].i=2;
			unit_cell.interaction[35].j=0;
			unit_cell.interaction[35].dx=1;
			unit_cell.interaction[35].dy=0;
			unit_cell.interaction[35].dz=1;
			//----------------------------------------------
			unit_cell.interaction[36].i=3;
			unit_cell.interaction[36].j=1;
			unit_cell.interaction[36].dx=-1;
			unit_cell.interaction[36].dy=0;
			unit_cell.interaction[36].dz=0;
			//----------------------------------------------
			unit_cell.interaction[37].i=3;
			unit_cell.interaction[37].j=2;
			unit_cell.interaction[37].dx=-1;
			unit_cell.interaction[37].dy=0;
			unit_cell.interaction[37].dz=0;
			//----------------------------------------------
			unit_cell.interaction[38].i=3;
			unit_cell.interaction[38].j=1;
			unit_cell.interaction[38].dx=-1;
			unit_cell.interaction[38].dy=0;
			unit_cell.interaction[38].dz=1;
			//----------------------------------------------
			unit_cell.interaction[39].i=3;
			unit_cell.interaction[39].j=2;
			unit_cell.interaction[39].dx=-1;
			unit_cell.interaction[39].dy=1;
			unit_cell.interaction[39].dz=0;
			//----------------------------------------------
			unit_cell.interaction[40].i=3;
			unit_cell.interaction[40].j=0;
			unit_cell.interaction[40].dx=0;
			unit_cell.interaction[40].dy=0;
			unit_cell.interaction[40].dz=0;
			//----------------------------------------------
			unit_cell.interaction[41].i=3;
			unit_cell.interaction[41].j=1;
			unit_cell.interaction[41].dx=0;
			unit_cell.interaction[41].dy=0;
			unit_cell.interaction[41].dz=0;
			//----------------------------------------------
			unit_cell.interaction[42].i=3;
			unit_cell.interaction[42].j=2;
			unit_cell.interaction[42].dx=0;
			unit_cell.interaction[42].dy=0;
			unit_cell.interaction[42].dz=0;
			//----------------------------------------------
			unit_cell.interaction[43].i=3;
			unit_cell.interaction[43].j=0;
			unit_cell.interaction[43].dx=0;
			unit_cell.interaction[43].dy=0;
			unit_cell.interaction[43].dz=1;
			//----------------------------------------------
			unit_cell.interaction[44].i=3;
			unit_cell.interaction[44].j=1;
			unit_cell.interaction[44].dx=0;
			unit_cell.interaction[44].dy=0;
			unit_cell.interaction[44].dz=1;
			//----------------------------------------------
			unit_cell.interaction[45].i=3;
			unit_cell.interaction[45].j=0;
			unit_cell.interaction[45].dx=0;
			unit_cell.interaction[45].dy=1;
			unit_cell.interaction[45].dz=0;
			//----------------------------------------------
			unit_cell.interaction[46].i=3;
			unit_cell.interaction[46].j=2;
			unit_cell.interaction[46].dx=0;
			unit_cell.interaction[46].dy=1;
			unit_cell.interaction[46].dz=0;
			//----------------------------------------------
			unit_cell.interaction[47].i=3;
			unit_cell.interaction[47].j=0;
			unit_cell.interaction[47].dx=0;
			unit_cell.interaction[47].dy=1;
			unit_cell.interaction[47].dz=1;
			//----------------------------------------------

		}
		else{
			std::cerr << "Error -  unknown crystal_type" << std::endl; 
			err::vexit();
		}

	}

} // end of namespace cs

