//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>

// Vampire headers
#include "create.hpp"
#include "errors.hpp"
#include "grains.hpp"
#include "material.hpp"
#include "random.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

// grains module headers
#include "internal.hpp"

namespace grains{

//====================================================================================
//
//														voronoi substructure
//
//    Function to generate a granular 2D substructure for a defined shape,
//    combining a voronoi construction with a bubble-like domain
//
//		(c) R F L Evans 17/11/2016
//
//             ______      _______
//           /       \    /       \
//          |        |   |        |
//          |        |   |        |
//       ---------------------------------------
//          |        |   |        |
//          |        |   |        |
//          \       /    \       /
//           ______       _______
//
//
//
//====================================================================================
//
void voronoi_substructure(std::vector<cs::catom_t> & catom_array){

	// No-op unless create:voronoi-grain-substructure was set
	if(!internal::generate_voronoi_substructure) return;

	//---------------------------------------------------
	// Local constants
	//---------------------------------------------------
	double grain_sd=internal::voronoi_sd;

	// Set number of particles in x and y directions
	double size = internal::voronoi_grain_substructure_size + internal::voronoi_grain_substructure_spacing;
	double grain_cell_size_x = size;
	double grain_cell_size_y = sqrt(3.0)*size;

	int num_x_particle = 4+vmath::iround(cs::system_dimensions[0]/(grain_cell_size_x));
	int num_y_particle = 4+vmath::iround(cs::system_dimensions[1]/(grain_cell_size_y));

	int init_num_grains = num_x_particle*num_y_particle*2;

	// Define initial grain arrays
	std::vector <std::vector <double> > grain_coord_array;
	std::vector <std::vector <std::vector <double> > > grain_vertices_array;

	// Reserve space for pointers
	grain_coord_array.reserve(init_num_grains);
	grain_vertices_array.reserve(init_num_grains);

	// Calculate pointers
	for(int grain=0;grain<init_num_grains;grain++){
		grain_coord_array.push_back(std::vector <double>());
		grain_coord_array[grain].reserve(2);
		grain_vertices_array.push_back(std::vector <std::vector <double> >());
		//for(int vertex=0;vertex<max_vertices;vertex++){
		//	grain_vertices_array[grain].push_back(std::vector <double>());
		//	grain_vertices_array[grain][vertex].reserve(2);
		//}
		//std::cout << grain_vertices_array[grain].size() << " " << grain_vertices_array[grain].capacity() << std::endl;
	}
	//std::cin.get();
	double delta_particle_x = grain_cell_size_x;
	double delta_particle_y = grain_cell_size_y;
	double delta_particle_x_parity = delta_particle_x*0.5;
	double delta_particle_y_parity = delta_particle_y*0.5;

	// Seed the grains module's own random generator so substructure
	// generation is reproducible independently of the global random stream.
	internal::grnd.seed(internal::grain_structure_seed);

	// Loop to generate hexagonal lattice points
	double particle_coords[2];

	int vp=int(internal::parity);
	int grain=0;

	for (int x_particle=0;x_particle < num_x_particle;x_particle++){
		for (int y_particle=0;y_particle < num_y_particle;y_particle++){
			for (int particle_parity=0;particle_parity<2;particle_parity++){

				//particle_coords[0] = (particle_parity)*delta_particle_x_parity + delta_particle_x*x_particle-size + vp*double(1-2*particle_parity)*delta_particle_x_parity;
				//particle_coords[1] = (particle_parity)*delta_particle_y_parity + delta_particle_y*y_particle-size;
				particle_coords[0] = (particle_parity)*delta_particle_x_parity + delta_particle_x*x_particle + vp*double(1-2*particle_parity)*delta_particle_x_parity;
				particle_coords[1] = (particle_parity)*delta_particle_y_parity + delta_particle_y*y_particle;

				grain_coord_array[grain].push_back(particle_coords[0]+grain_sd*mtrandom::gaussianc(internal::grnd)*delta_particle_x);
				grain_coord_array[grain].push_back(particle_coords[1]+grain_sd*mtrandom::gaussianc(internal::grnd)*delta_particle_y);

				grain++;
			}
		}
	}
	//-----------------------
	// Check for grains >=1
	//-----------------------
	if(grain<1){
		terminaltextcolor(RED);
		std::cerr << "Error! - No grains found in structure - Increase system dimensions" << std::endl;
		terminaltextcolor(WHITE);
		zlog << zTs() << "Error! - No grains found in structure - Increase system dimensions" << std::endl;
		err::vexit();
	}

	// Tessellate the hexagonal lattice of substructure sites into an ordinary
	// (uniform-weight) Voronoi diagram, clipped to the film's footprint.
	internal::populate_vertex_points_power(grain_coord_array, grain_vertices_array, std::vector<double>(),
		cs::system_dimensions[0], cs::system_dimensions[1], false, false);

	// Pull each cell boundary inward by half the requested substructure
	// spacing, giving neighbouring sub-grains a constant absolute gap
	// regardless of their individual cell size.
	internal::apply_grain_spacing(grain_coord_array, grain_vertices_array,
		0.5*internal::voronoi_grain_substructure_spacing);

   // round grain corners if requested (create:grain-rounding; 0.0 = off)
	if(internal::grain_rounding > 0.0){
		internal::voronoi_grain_rounding(grain_coord_array, grain_vertices_array, internal::grain_rounding);
	}

   // Bin atoms into supercells once, shared by every grain's atom-assignment pass
   const internal::supercell_array_t supercell_array =
      internal::build_supercell_array(catom_array, cs::unit_cell.dimensions[0], cs::unit_cell.dimensions[1],
                                               cs::total_num_unit_cells[0], cs::total_num_unit_cells[1]);

	std::cout <<"Generating voronoi substructure";
	zlog << zTs() << "Generating voronoi substructure";

   // array to store if atoms are included in substructure (assume not)
   std::vector<bool> insub(catom_array.size(),false);

   const double ssz = cs::system_dimensions[2];

   //------------------------------------------------------------------------------
   // Give each 2D substructure cell a height-dependent in-plane radius, so
   // each sub-grain bulges out from a nucleation height and tapers back in
   // toward the top and bottom surfaces - a lens/ellipsoidal cross-section
   // rather than a straight-sided column. factor_radius below scales the
   // in-plane point-in-polygon test radius as a function of the reduced
   // height frh=z/thickness relative to the nucleation height, with exponent
   // sphere_radius controlling how sharply the radius falls off.
   //------------------------------------------------------------------------------
   const double sphere_radius = internal::voronoi_grain_substructure_crystallization_radius; //*internal::voronoi_grain_substructure_size*radius_factor;

	// loop over all grains with vertices
	for(unsigned int grain=0;grain<grain_coord_array.size();grain++){
		// Exclude grains with zero vertices
		internal::print_grain_progress(grain, grain_coord_array.size());
		if(grain_vertices_array[grain].size()!=0){

         // determine coordinate offset for grains
         const double x0 = grain_coord_array[grain][0];
         const double y0 = grain_coord_array[grain][1];

         internal::grain_footprint_t fp =
            internal::compute_grain_footprint(grain_vertices_array[grain], x0, y0,
                                                       cs::unit_cell.dimensions[0], cs::unit_cell.dimensions[1]);
         const int num_vertices = fp.px.size();

         // copy overlap of substructure grains to a local constant
         const double overlap = internal::voronoi_grain_substructure_overlap_factor;

         internal::assign_atoms_in_footprint(supercell_array, fp, [&](int atom){

            // Get atomic position
            const double x = catom_array[atom].x;
            const double y = catom_array[atom].y;
            const double z = catom_array[atom].z;
            const double frh = z/ssz;
            const double nucleation_height = create::get_material_substructure_nucleation_height(catom_array[atom].material);

            int mat = catom_array[atom].material;
            // calculate reduced ranges for materials with small offset to prevent dangling atoms
            double rminz = create::get_material_height_min(mat)-0.01;
            double rmaxz = create::get_material_height_max(mat)+0.01;
            double factor_radius = 0.0;
            if(frh > nucleation_height){
               // multiply by small factor to ensure grains touch at boundary for zero spacing
               factor_radius = 1.04*pow(1.0+((nucleation_height-frh)/(rmaxz-nucleation_height)),sphere_radius);
            }
            else{
               factor_radius = 1.04*pow((1.0-(frh-nucleation_height)/(rminz-nucleation_height)),sphere_radius);
            }
            if(vmath::point_in_polygon_scaled(x-x0,y-y0,overlap*factor_radius,fp.px.data(),fp.py.data(),num_vertices)){
               insub[atom] = true;
            }

         });
		}
	}

	terminaltextcolor(GREEN);
	std::cout << "done!" << std::endl;
	terminaltextcolor(WHITE);
	zlog << "done!" << std::endl;

   // Now fill in with fill materials
   for(int mat=0;mat<mp::num_materials;mat++){
      if(create::get_material_sub_fill(mat)){
         double min = create::get_material_height_min(mat)*cs::system_dimensions[2];
         double max = create::get_material_height_max(mat)*cs::system_dimensions[2];

         // loop over all atoms selecting only deselected atoms within min/max
         for(unsigned int atom=0;atom<catom_array.size();atom++){
            if( (catom_array[atom].z < max) && (catom_array[atom].z >= min) && (catom_array[atom].include==true && insub[atom] == false)){
               // set atom to fill material
               catom_array[atom].material=mat;
               // re-include atom
               insub[atom] = true;
            }
         }
      }
   }


   // Now delete atoms not in substructure
   for(unsigned int atom=0; atom < catom_array.size(); atom++){
      if(insub[atom] == false) catom_array[atom].include=false;
   }

	// check for continuous layer
	for(unsigned int atom=0; atom < catom_array.size(); atom++){
	  if(mp::material[catom_array[atom].material].continuous==true && catom_array[atom].include == false ){
	    catom_array[atom].include=true;
	    catom_array[atom].grain=int(grain_coord_array.size()-1);
	  }
	}

	return;
}

} // end of namespace grains
