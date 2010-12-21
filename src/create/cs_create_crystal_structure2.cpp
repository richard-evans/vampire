//======================================================================
//                         create_crystal_structure
//   Subroutine to set system size and create desired crystal structure
//
//======================================================================

#include "errors.hpp"
#include "create.hpp"
#include "material.hpp"
#include "vmpi.hpp"
#include <string>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>

namespace cs{
	
int create_crystal_structure(std::vector<cs::catom_t> & catom_array){
	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "cs::create_crystal_structure has been called" << std::endl;}	

	//-----------------------------------------------------------------------
	// Pre Determine crystal structure and estimate number of atoms each material
	//-----------------------------------------------------------------------
	std::vector<int> mat_min(mp::num_materials);
	std::vector<int> mat_max(mp::num_materials);
	std::vector<std::string> mat_cs(mp::num_materials);
	
	const int max_atoms_per_unit_cell=4;
	int supercell_dim[3];

	supercell_dim[0] = iround(mp::system_dimensions[0]/mp::lattice_constant[0]);
	supercell_dim[1] = iround(mp::system_dimensions[1]/mp::lattice_constant[1]);
	supercell_dim[2] = iround(mp::system_dimensions[2]/mp::lattice_constant[2]);
	
	int ncell[3];
	
	#ifdef MPICF
		ncell[0] = iround((vmpi::max_dimensions[0]-vmpi::min_dimensions[0])/mp::lattice_constant[0]);
		ncell[1] = iround((vmpi::max_dimensions[1]-vmpi::min_dimensions[1])/mp::lattice_constant[1]);
		ncell[2] = iround((vmpi::max_dimensions[2]-vmpi::min_dimensions[2])/mp::lattice_constant[2]);
	#else
		ncell[0]=supercell_dim[0];
		ncell[1]=supercell_dim[1];
		ncell[2]=supercell_dim[2];
	#endif
	
	int num_atoms=ncell[0]*ncell[1]*ncell[2]*max_atoms_per_unit_cell;

	for(int mat=0;mat<mp::num_materials;mat++){
		mat_cs[mat] = mp::material[mat].crystal_structure;
		// supercell z min and max
		mat_min[mat]=iround(mp::material[mat].min*mp::system_dimensions[2]/mp::lattice_constant[2]);
		mat_max[mat]=int(mp::material[mat].max*mp::system_dimensions[2]/mp::lattice_constant[2])+1; //always round up number of cells
		//std::cout << mat_min[mat] << "\t" << mat_max[mat] << std::endl;
		
	}

	// set catom_array size
	catom_array.reserve(num_atoms);
	// Initialise atoms number
	int atom=0;

	//-----------------------------------------------------------------------
	// Create crystal structure for each material
	//-----------------------------------------------------------------------
	for(int mat=0;mat<mp::num_materials;mat++){
		
		//--------------------------------------------------------------
		// Determine atoms per unit cell and unit cell atom coordinates
		//--------------------------------------------------------------
		int atoms_per_unit_cell=0;
		std::vector <double> unit_cell_coords(0);
		
		// Simple Cubic
		if(mat_cs[mat]=="sc"){
			atoms_per_unit_cell=1;
			unit_cell_coords.resize(3*atoms_per_unit_cell);
			//-----------------------------
			unit_cell_coords[0]=0.0;
			unit_cell_coords[1]=0.0;
			unit_cell_coords[2]=0.0;
			//-----------------------------
		}
		if(mat_cs[mat]=="bcc"){
			atoms_per_unit_cell=2;
			unit_cell_coords.resize(3*atoms_per_unit_cell);
			//-----------------------------
			unit_cell_coords[0]=0.0;
			unit_cell_coords[1]=0.0;
			unit_cell_coords[2]=0.0;
			//-----------------------------
			unit_cell_coords[3]=0.5;
			unit_cell_coords[4]=0.5;
			unit_cell_coords[5]=0.5;
			//-----------------------------
		}
		if(mat_cs[mat]=="fct"){
			atoms_per_unit_cell=2;
			unit_cell_coords.resize(3*atoms_per_unit_cell);
			//-----------------------------
			unit_cell_coords[0]=0.0;
			unit_cell_coords[1]=0.0;
			unit_cell_coords[2]=0.0;
			//-----------------------------
			unit_cell_coords[3]=0.5;
			unit_cell_coords[4]=0.5;
			unit_cell_coords[5]=0.0;
			//-----------------------------
		}
		if(mat_cs[mat]=="fcc"){
			atoms_per_unit_cell=4;
			unit_cell_coords.resize(3*atoms_per_unit_cell);
			//-----------------------------
			unit_cell_coords[0]=0.0;
			unit_cell_coords[1]=0.0;
			unit_cell_coords[2]=0.0;
			//-----------------------------
			unit_cell_coords[3]=0.5;
			unit_cell_coords[4]=0.5;
			unit_cell_coords[5]=0.0;
			//-----------------------------
			unit_cell_coords[6]=0.5;
			unit_cell_coords[7]=0.0;
			unit_cell_coords[8]=0.5;
			//-----------------------------
			unit_cell_coords[9]=0.0;
			unit_cell_coords[10]=0.5;
			unit_cell_coords[11]=0.5;
			//-----------------------------
		}
		
		if(atoms_per_unit_cell==0){std::cout << "Error determining atoms_per_unit_cell, unknown crystal_type" << std::endl; exit(1);}

		double min=mp::material[mat].min*mp::system_dimensions[2];
		double max=mp::material[mat].max*mp::system_dimensions[2];
		
		int min_bounds[3];
		int max_bounds[3];
		
		#ifdef MPICF
			min_bounds[0] = int(vmpi::min_dimensions[0]/mp::lattice_constant[0]);
			min_bounds[1] = int(vmpi::min_dimensions[1]/mp::lattice_constant[1]);
			min_bounds[2] = int(vmpi::min_dimensions[2]/mp::lattice_constant[2]);
			max_bounds[0] = 1+int(vmpi::max_dimensions[0]/mp::lattice_constant[0]);
			max_bounds[1] = 1+int(vmpi::max_dimensions[1]/mp::lattice_constant[1]);
			max_bounds[2] = 1+int(vmpi::max_dimensions[2]/mp::lattice_constant[2]);
		#else
			min_bounds[0]=0;
			min_bounds[1]=0;
			min_bounds[2]=0;
			max_bounds[0]=supercell_dim[0];
			max_bounds[1]=supercell_dim[1];
			max_bounds[2]=supercell_dim[2];
		#endif
		
		// Duplicate unit cell
		for(int z=min_bounds[2];z<max_bounds[2];z++){
			for(int y=min_bounds[1];y<max_bounds[1];y++){
				for(int x=min_bounds[0];x<max_bounds[0];x++){
					// Loop over atoms in unit cell
					for(int uca=0;uca<atoms_per_unit_cell;uca++){
						double cx = (double(x)+unit_cell_coords[3*uca+0])*mp::lattice_constant[0];
						double cy = (double(y)+unit_cell_coords[3*uca+1])*mp::lattice_constant[1];
						double cz = (double(z)+unit_cell_coords[3*uca+2])*mp::lattice_constant[2];
						if((cz>=min) && (cz<max)){
							#ifdef MPICF
								if(	(cx>=vmpi::min_dimensions[0] && cx<vmpi::max_dimensions[0]) &&
										(cy>=vmpi::min_dimensions[1] && cy<vmpi::max_dimensions[1]) &&
										(cz>=vmpi::min_dimensions[2] && cz<vmpi::max_dimensions[2])){
							#endif
								catom_array.push_back(cs::catom_t());
								catom_array[atom].x=cx;
								catom_array[atom].y=cy;
								catom_array[atom].z=cz;
								catom_array[atom].material=mat;
								atom++;
							#ifdef MPICF
								}
							#endif
						}
					}
				}
			}
		}
	}
	// Check to see if actual nad expected number of atoms agree, if not trim the excess
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
		std::cout << "Error - no atoms have been generated, increase system dimensions!" << std::endl;
		exit(1);
	}

	return EXIT_SUCCESS;
}

}
