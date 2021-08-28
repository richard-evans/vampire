//-----------------------------------------------------------------------------
//
//  Vampire - A code for atomistic simulation of magnetic materials
//
//  Copyright (C) 2009-2012 R.F.L.Evans
//
//  Email:richard.evans@york.ac.uk
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
//  General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software Foundation,
//  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
//
// ----------------------------------------------------------------------------
//
//======================================================================
//                         create_system_type
//   Subroutine to set system size and create desired crystal structure
//
//======================================================================

// C++ standard library headers
#include <string>
#include <sstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <list>

// Vampire headers
#include "errors.hpp"
#include "create.hpp"
#include "grains.hpp"
#include "material.hpp"
#include "random.hpp"
#include "unitcell.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

// Internal create header
#include "internal.hpp"

namespace create{

   //----------------------------------------
   // function prototypes
   //----------------------------------------

   int intermixing(std::vector<cs::catom_t> &);
   void dilute(std::vector<cs::catom_t> &);
   void geometry(std::vector<cs::catom_t> &);
   void fill(std::vector<cs::catom_t> &);
   void roughness(std::vector<cs::catom_t> &);

//======================================================================
//                         create_system_type
//   Subroutine to set system size and create desired crystal structure
//
//======================================================================
int create_system_type(std::vector<cs::catom_t> & catom_array){

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "cs::create_system_type has been called" << std::endl;}

	//----------------------------------------
	// Local variables
	//----------------------------------------

	//----------------------------------------------------------------------------------
	// Choose which system type to create
	//----------------------------------------------------------------------------------
	switch(cs::system_creation_flags[2]){
		case 0: // Isolated particle
			internal::particle(catom_array);
			break;

		case 1: // Cubic Particle Array
			internal::particle_array(catom_array);
			break;

		case 2: // Hexagonal Particle Array
			internal::hex_particle_array(catom_array);
			break;

		case 3: // Voronoi Granular Film
			voronoi_film(catom_array);
			break;

		case 4: // Grain Growth Method
			//grain_growth(cs_num_atoms,cs_coord_array,particle_include_array,cs_atom_type_array);
			std::cerr << "Grain growth not yet implemented, exiting" << std::endl;
			err::vexit();
			break;

		default:{
			std::cerr << "Unknown system type requested, exiting" << std::endl;
			err::vexit();
		}
	}

   // Check for voronoi construction and apply before csg operations
   if(create::internal::generate_voronoi_substructure) create::internal::voronoi_substructure(catom_array);

	// call fill function to fill in void
	fill(catom_array);

	// call geometry function
	geometry(catom_array);

	// call intermixing function - must be before alloy function
	intermixing(catom_array);

	// call surface roughness function
	// Formally called here but now moved to crystal structure generation
	//roughness(catom_array);

	// call alloy function
	create::internal::alloy(catom_array);

	// call dilution function
	dilute(catom_array);

	// Delete unneeded atoms
	create::internal::clear_atoms(catom_array);

	// Calculate final atomic composition
	create::internal::calculate_atomic_composition(catom_array);

   // For parallel check which processors have zero atoms
   #ifdef MPICF
      uint64_t num_atoms_check = 0;
      // if a processor has zero atoms then flag as 1 (0 has more than zero atoms)
      if(catom_array.size() == 0 ) num_atoms_check = 1;
      // Check globally for no errors
      MPI_Allreduce(MPI_IN_PLACE, &num_atoms_check, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
      // If error, determine which ranks have no atoms
      if( num_atoms_check > 0){
         std::vector<uint64_t> no_atoms(vmpi::num_processors, 0);
         if(catom_array.size() == 0 ) no_atoms[vmpi::my_rank] = 1;
         MPI_Allreduce( MPI_IN_PLACE , &no_atoms[0], vmpi::num_processors, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
         // generate error message
         std::stringstream message_stream;
         if(vmpi::my_rank == 0){
            message_stream << "Error! the following processes have zero atoms: ";
            for(int p=0; p < vmpi::num_processors; p++){
               if(no_atoms[p] == 1) message_stream << p << " ";
            }
            message_stream << ". All parallel processes must contain atoms - change system dimensions or add fill material!";
         }
         // Output error message to screen
         terminaltextcolor(RED);
            std::cout << "Error, no atoms generated on some processors for requested system shape - change system dimensions or add fill material!" << std::endl;
         terminaltextcolor(WHITE);
         // Exit without abort and nice error message
         err::v_parallel_all_exit(message_stream.str());
      }
   #else
		// Check for zero atoms generated
		if(catom_array.size()==0){
			terminaltextcolor(RED);
			std::cerr << "Error, no atoms generated for requested system shape - increase system dimensions or reduce particle size!" << std::endl;
			terminaltextcolor(WHITE);
			err::vexit();
		}
   #endif

	return 0;
}

/*
int pop_template_2D(int** template_array_2D,int cs_num_atoms,int** cs_coord_array,int* particle_include_array,int* template_bounds){

	//====================================================================================
	//
	//											pop_template_2D
	//
	//					Subroutine to map 2D template onto particles include array
	//
	//									Version 0.2 R Evans 16/09/2008
	//
	//====================================================================================

  	for(int atom=0;atom<cs_num_atoms;atom++){
		#ifdef MPICF
		int cx=cs_coord_array[atom][0]+mpi_create_variables::int_mpi_offset[0];
		int cy=cs_coord_array[atom][1]+mpi_create_variables::int_mpi_offset[1];
		#else
		int cx=cs_coord_array[atom][0];
		int cy=cs_coord_array[atom][1];
		#endif

		// Check for coords in template bounds
		if((cx>0) && (cx<template_bounds[0]) && (cy>0) && (cy<template_bounds[1])){
			if(template_array_2D[cx][cy]!=-1){
				particle_include_array[atom]=template_array_2D[cx][cy];
			}
		}
	}

	return 0;

}
*/
namespace internal{

void clear_atoms(std::vector<cs::catom_t> & catom_array){

   // check calling of routine if error checking is activated
   if(err::check==true){std::cout << "cs::clear_atoms has been called" << std::endl;}

   // Get original and new number of atoms
   const int num_atoms=catom_array.size();
   int num_included=0;
   for(int a=0;a<num_atoms;a++){
      if(catom_array[a].include == true && mp::material[catom_array[a].material].non_magnetic != 1){
         num_included++;
      }
   }

   // Get total number of atoms contributing to the body volume, i.e. all the atoms generated within the body shape
   create::num_total_atoms_non_filler = 0;
   for(int a=0;a<num_atoms;a++){
      if(catom_array[a].include == true && mp::material[catom_array[a].material].fill == false && mp::material[catom_array[a].material].non_magnetic == 1){
         create::num_total_atoms_non_filler++;
      }
   }

   // check if there are unneeded atoms
   if(num_atoms!=num_included){
      // create temporary copy for atoms
      std::vector<cs::catom_t> tmp_catom_array(num_atoms);
      tmp_catom_array=catom_array;
      // resize original array to new number of atoms
      catom_array.resize(num_included);
      int atom=0;
      // loop over all existing atoms
      for(int a=0;a<num_atoms;a++){
         // if atom is to be included and is non-magnetic copy to new array
         if(catom_array[a].include==true && mp::material[catom_array[a].material].non_magnetic != 1 ){
            catom_array[atom]=tmp_catom_array[a];
            atom++;
         }
         // if atom is part of a non-magnetic material to be removed then save to nm array
         else if(catom_array[a].include == true && mp::material[catom_array[a].material].non_magnetic == 1){
            cs::nm_atom_t tmp;
         	tmp.x = catom_array[a].x;
         	tmp.y = catom_array[a].y;
         	tmp.z = catom_array[a].z;
         	tmp.mat = catom_array[a].material;
            tmp.cat = catom_array[atom].lh_category;
         	// save atom to non-magnet array
         	cs::non_magnetic_atoms_array.push_back(tmp);
         }
      }
      tmp_catom_array.resize(0);

      zlog << zTs() << "Removed " << cs::non_magnetic_atoms_array.size() << " non-magnetic atoms from system" << std::endl;

   }

   return;

}
}

int intermixing(std::vector<cs::catom_t> & catom_array){
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "cs::intermixing has been called" << std::endl;}

   // re-seed random number generator on each CPU with a different number
	create::internal::grnd.seed(vmpi::parallel_rng_seed(create::internal::mixing_seed));

	// loop over all atoms
	for(unsigned int atom=0;atom<catom_array.size();atom++){
		// get current material
		int current_material=catom_array[atom].material;
		int final_material=current_material;

		//loop over all potential intermixing materials
		for(int mat=0;mat<mp::num_materials;mat++){
			if(mp::material[current_material].intermixing[mat]>0.0){
				// find which region atom is in and test for probability of different material
				double z=catom_array[atom].z;
				double min = create::internal::mp[mat].min*cs::system_dimensions[2];
				double max = create::internal::mp[mat].max*cs::system_dimensions[2];
				double mean = (min+max)/2.0;
				if(z<=min){
					double probability=0.5+0.5*tanh((z-min)/(mp::material[current_material].intermixing[mat]*cs::system_dimensions[2]));
					if(create::internal::grnd() < probability) final_material=mat;
				}
				else if(z>min && z<=mean){
					double probability=0.5+0.5*tanh((z-min)/(mp::material[current_material].intermixing[mat]*cs::system_dimensions[2]));
					if(create::internal::grnd() < probability) final_material=mat;
				}
				else if(z>mean && z<=max){
					double probability=0.5-0.5*tanh((z-max)/(mp::material[current_material].intermixing[mat]*cs::system_dimensions[2]));
					if(create::internal::grnd() < probability) final_material=mat;
				}
				else if(z>max){
					double probability=0.5-0.5*tanh((z-max)/(mp::material[current_material].intermixing[mat]*cs::system_dimensions[2]));
					//std::cout << current_material << "\t" << mat << "\t" << atom << "\t" << z << "\t" << max << "\t" << probability << std::endl;
					if(create::internal::grnd() < probability) final_material=mat;
				}
			}
		}

		// set final material
		catom_array[atom].material=final_material;
	}

	return EXIT_SUCCESS;
}

void dilute (std::vector<cs::catom_t> & catom_array){
   // check calling of routine if error checking is activated
   if(err::check==true){std::cout << "cs::dilute has been called" << std::endl;}

   // re-seed random number generator on each CPU with a different number
   create::internal::grnd.seed(vmpi::parallel_rng_seed(create::internal::dilute_seed));

   // loop over all atoms
   for(unsigned int atom=0;atom<catom_array.size();atom++){
      // if atom material is alloy master
      int local_material=catom_array[atom].material;
      double probability = mp::material[local_material].density;
      if(create::internal::grnd() > probability) catom_array[atom].include=false;
   }

   return;
}

void geometry (std::vector<cs::catom_t> & catom_array){
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "cs::geometry has been called" << std::endl;}

	// Check for any geometry
	bool cut=false;

	for(int mat=0; mat<mp::num_materials; mat++){
		if(mp::material[mat].geometry>0) cut=true;
	}

	// Return from function if no geometry is defined.
	if(cut==false) return;

	// Otherwise proceed
	zlog << zTs() << "Cutting materials within defined geometry." << std::endl;

	// Check for force material type by geometry
	if(cs::SelectMaterialByGeometry==false){

		// loop over all atoms
		for(unsigned int atom=0;atom<catom_array.size();atom++){

			// check for geometry information
			const int geo=mp::material[catom_array[atom].material].geometry;

			// if exists, then remove atoms outside polygon
			if(geo!=0){
				double x = catom_array[atom].x;
				double y = catom_array[atom].y;
				std::vector<double> px(geo);
				std::vector<double> py(geo);
				// Initialise polygon points
				for(int p=0;p<geo;p++){
					px[p]=mp::material[catom_array[atom].material].geometry_coords[p][0]*cs::system_dimensions[0];
					py[p]=mp::material[catom_array[atom].material].geometry_coords[p][1]*cs::system_dimensions[1];
				}
				// check if point is outside of polygon, if so delete it
				if(vmath::point_in_polygon2(x,y,px,py,geo)==false){
					catom_array[atom].include=false;
				}
			}

		} 

	}
	else{

		// Re-identify all atoms as material 0 and exclude by default
		for(unsigned int atom=0;atom<catom_array.size();atom++){
			catom_array[atom].material=0;
			catom_array[atom].include=false;
		}

      std::cout << "AAAA" <<std::endl;

		// loop over all materials and include according to geometry

		// determine z-bounds for materials
		std::vector<double> mat_min(mp::num_materials);
		std::vector<double> mat_max(mp::num_materials);
    std::vector<double> mat_min_x(mp::num_materials);
		std::vector<double> mat_max_x(mp::num_materials);

		// initialise array to hold z-bound information in real space
		for(int mat=0;mat<mp::num_materials;mat++){
			mat_min[mat]=create::internal::mp[mat].min*cs::system_dimensions[2];
			mat_max[mat]=create::internal::mp[mat].max*cs::system_dimensions[2];
      mat_min_x[mat]=0.0*cs::system_dimensions[0];
      mat_max_x[mat]=0.5*cs::system_dimensions[0];
			// alloys generally are not defined by height, and so have max = 0.0
			if(mat_max[mat]<0.0000001) mat_max[mat]=-0.1;
		}

		// Loop over all materials
		for(int mat=0;mat<mp::num_materials;mat++){

			// check for geometry information
			const int geo=mp::material[mat].geometry;

			// if geometry information exists, then include atoms inside polygon and within material height
			if(geo!=0){

				// create array to store geoemtric points
				std::vector<double> px(geo);
				std::vector<double> py(geo);
				// Initialise polygon points
				for(int p=0;p<geo;p++){
					px[p]=mp::material[mat].geometry_coords[p][0]*cs::system_dimensions[0];
					py[p]=mp::material[mat].geometry_coords[p][1]*cs::system_dimensions[1];
				}

				for(unsigned int atom=0;atom<catom_array.size();atom++){
					double x = catom_array[atom].x;
					double y = catom_array[atom].y;
					const double z = catom_array[atom].z;
          std::cout << x << '\t' << mat_min_x[mat] << '\t' << mat_max_x[mat] <<std::endl;
					if((z>=mat_min[mat]) && (z<mat_max[mat]) && (x>=mat_min_x[mat]) && (x<mat_max_x[mat]) && (vmath::point_in_polygon2(x,y,px,py,geo)==true)){
						catom_array[atom].material=mat;
						catom_array[atom].include=true;
					}
				}
			}
		}
	}

	return;
}

//-----------------------------------------------------------
//
///  Function to replace deleted atoms with in-fill material
//
///  Can be used to create embedded nanoparticle arrays or
///  granular recording media.
//
///  v1 18/09/2013
///  (c) R F L Evans
//
//-----------------------------------------------------------
void fill(std::vector<cs::catom_t> & catom_array){
   // check calling of routine if error checking is activated
   if(err::check==true){std::cout << "cs::fill has been called" << std::endl;}

   //loop over all potential fill materials
   for(int mat=0;mat<mp::num_materials;mat++){

      // If material is fill material
      if(mp::material[mat].fill){

         // Calculate min/max heights for fill material (for a partial fill)
         double min = create::internal::mp[mat].min*cs::system_dimensions[2];
         double max = create::internal::mp[mat].max*cs::system_dimensions[2];

         // loop over all atoms selecting only deselected atoms within min/max
         for(unsigned int atom=0;atom<catom_array.size();atom++){
            if( (catom_array[atom].z < max) && (catom_array[atom].z >= min) && (catom_array[atom].include==false)){
               // set atom to fill material
               catom_array[atom].material=mat;
               // re-include atom
               catom_array[atom].include=true;
            }
         }

      }

   }

   return;

}



} // end of namespace create
