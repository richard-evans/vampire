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
#include "create.hpp"
#include "errors.hpp"
#include "grains.hpp"
#include "material.hpp"
#include "qvoronoi.hpp"
#include "random.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "vmath.hpp"
#include "voronoi.hpp"

#include <cmath>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>

#include "internal.hpp"

namespace create_voronoi{
	bool parity=0;	/// left-right (0) or right-left (1) point initialisation
	bool rounded=false;
	double area_cutoff=0.8;
	double voronoi_sd=0.15;			/// Standard Deviation of voronoi grains

}

namespace cs{

int voronoi_film(std::vector<cs::catom_t> & catom_array){

	// check calling of routine if error checking is activated
	if(err::check==true){
      terminaltextcolor(RED);
      std::cerr << "cs::voronoi_film has been called" << std::endl;
      terminaltextcolor(WHITE);
   }
	//====================================================================================
	//
	//														voronoi
	//
	//				Subroutine to create granular system using qhull voronoi generator
	//
	//							Version 1.0 R Evans 16/07/2009
	//
	//====================================================================================
	//
	//		Locally allocated variables: 	init_grain_coord_array
	//												init_grain_pointx_array
	//												init_grain_pointy_array
	//												init_num_assoc_vertices_array
	//
	//=====================================================================================


	//---------------------------------------------------
	// Local constants
	//---------------------------------------------------
	const int max_vertices=50;
	double grain_sd=create_voronoi::voronoi_sd;

	// Set number of particles in x and y directions
	double size = create::internal::voronoi_grain_size + create::internal::voronoi_grain_spacing;
	double grain_cell_size_x = size;
	double grain_cell_size_y = sqrt(3.0)*size;

	int num_x_particle = 4+2*vmath::iround(cs::system_dimensions[0]/(grain_cell_size_x));
	int num_y_particle = 4+2*vmath::iround(cs::system_dimensions[1]/(grain_cell_size_y));

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

	// Set voronoi seed;
	mtrandom::grnd.seed(mtrandom::voronoi_seed);

	// Loop to generate hexagonal lattice points
	double particle_coords[2];

	int vp=int(create_voronoi::parity);
	int grain=0;

	for (int x_particle=0;x_particle < num_x_particle;x_particle++){
		for (int y_particle=0;y_particle < num_y_particle;y_particle++){
			for (int particle_parity=0;particle_parity<2;particle_parity++){

				//particle_coords[0] = (particle_parity)*delta_particle_x_parity + delta_particle_x*x_particle-size + vp*double(1-2*particle_parity)*delta_particle_x_parity;
				//particle_coords[1] = (particle_parity)*delta_particle_y_parity + delta_particle_y*y_particle-size;
				particle_coords[0] = (particle_parity)*delta_particle_x_parity + delta_particle_x*x_particle + vp*double(1-2*particle_parity)*delta_particle_x_parity;
				particle_coords[1] = (particle_parity)*delta_particle_y_parity + delta_particle_y*y_particle;

				grain_coord_array[grain].push_back(particle_coords[0]+grain_sd*mtrandom::gaussian()*delta_particle_x);
				grain_coord_array[grain].push_back(particle_coords[1]+grain_sd*mtrandom::gaussian()*delta_particle_y);

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

   // Calculate Voronoi construction using qhull removing boundary grains (false)
	create::internal::populate_vertex_points(grain_coord_array, grain_vertices_array, false);

	// Shrink Voronoi vertices in reduced coordinates to get spacing
	double shrink_factor = create::internal::voronoi_grain_size/(create::internal::voronoi_grain_size+create::internal::voronoi_grain_spacing);

	// Reduce vertices to relative coordinates
	for(unsigned int grain=0;grain<grain_coord_array.size();grain++){
		//grain_coord_array[grain][0]=0.0;
		//grain_coord_array[grain][1]=0.0;
		const int nv = grain_vertices_array[grain].size();
		// Exclude grains with zero vertices
      if(nv!=0){
			for(int vertex=0;vertex<nv;vertex++){
				double vc[2]; // vertex coordinates
				vc[0]=shrink_factor*(grain_vertices_array[grain][vertex][0]-grain_coord_array[grain][0]);
				vc[1]=shrink_factor*(grain_vertices_array[grain][vertex][1]-grain_coord_array[grain][1]);
				grain_vertices_array[grain][vertex][0]=vc[0];
				grain_vertices_array[grain][vertex][1]=vc[1];
			}
		}
	}

   // round grains if necessary
	if(create_voronoi::rounded) create::internal::voronoi_grain_rounding(grain_coord_array, grain_vertices_array);

	// Create a 2D supercell array of atom numbers to improve performance for systems with many grains
	std::vector < std::vector < std::vector < int > > > supercell_array;

	int min_bounds[3];
	int max_bounds[3];

	min_bounds[0]=0;
	min_bounds[1]=0;
	min_bounds[2]=0;
	max_bounds[0]=cs::total_num_unit_cells[0];
	max_bounds[1]=cs::total_num_unit_cells[1];
	max_bounds[2]=cs::total_num_unit_cells[2];

	// allocate supercell array
	int dx = max_bounds[0]-min_bounds[0];
	int dy = max_bounds[1]-min_bounds[1];

	supercell_array.resize(dx);
	for(int i=0;i<dx;i++) supercell_array[i].resize(dy);

	// loop over atoms and populate supercell array
	for(unsigned int atom=0;atom<catom_array.size();atom++){
		int cx = int (catom_array[atom].x/unit_cell.dimensions[0]);
		int cy = int (catom_array[atom].y/unit_cell.dimensions[1]);
		supercell_array.at(cx).at(cy).push_back(atom);
	}

	// Determine order for core-shell grains
   std::list<create::internal::core_radius_t> material_order(0);
   for(int mat=0;mat<mp::num_materials;mat++){
      create::internal::core_radius_t tmp;
      tmp.mat=mat;
      tmp.radius=mp::material[mat].core_shell_size;
      material_order.push_back(tmp);
   }
   // sort by increasing radius
   material_order.sort(create::internal::compare_radius);

	std::cout <<"Generating Voronoi Grains";
	zlog << zTs() << "Generating Voronoi Grains";

   // arrays to store list of grain vertices
   double tmp_grain_pointx_array[max_vertices];
	double tmp_grain_pointy_array[max_vertices];

	// loop over all grains with vertices
	for(unsigned int grain=0;grain<grain_coord_array.size();grain++){
		// Exclude grains with zero vertices

		if((grain%(grain_coord_array.size()/10))==0){
		  std::cout << "." << std::flush;
		  zlog << "." << std::flush;
		}
		if(grain_vertices_array[grain].size()!=0){

			// initialise minimum and max supercell coordinates for grain
			int minx=10000000;
			int maxx=0;
			int miny=10000000;
			int maxy=0;

			// Set temporary vertex coordinates (real) and compute cell ranges
			int num_vertices = grain_vertices_array[grain].size();
			for(int vertex=0;vertex<num_vertices;vertex++){
				// determine vertex coordinates
				tmp_grain_pointx_array[vertex]=grain_vertices_array[grain][vertex][0];
				tmp_grain_pointy_array[vertex]=grain_vertices_array[grain][vertex][1];
				// determine unit cell coordinates encompassed by grain
				int x = int((tmp_grain_pointx_array[vertex]+grain_coord_array[grain][0])/unit_cell.dimensions[0]);
				int y = int((tmp_grain_pointy_array[vertex]+grain_coord_array[grain][1])/unit_cell.dimensions[1]);
				if(x < minx) minx = x;
				if(x > maxx) maxx = x;
				if(y < miny) miny = y;
				if(y > maxy) maxy = y;
			}

			// determine coordinate offset for grains
			const double x0 = grain_coord_array[grain][0];
			const double y0 = grain_coord_array[grain][1];

			// loopover cells
			for(int i=minx;i<=maxx;i++){
				for(int j=miny;j<=maxy;j++){

					// loop over atoms in cells;
					for(unsigned int id=0;id<supercell_array[i][j].size();id++){
						int atom = supercell_array[i][j][id];

						// Get atomic position
						double x = catom_array[atom].x;
						double y = catom_array[atom].y;

						if(mp::material[catom_array[atom].material].core_shell_size>0.0){
							// Iterate over materials
							for(std::list<create::internal::core_radius_t>::iterator it = material_order.begin(); it !=  material_order.end(); it++){
								int mat = (it)->mat;
								double factor = mp::material[mat].core_shell_size;
								double maxz=create::internal::mp[mat].max*cs::system_dimensions[2];
								double minz=create::internal::mp[mat].min*cs::system_dimensions[2];
								double cz=catom_array[atom].z;
                        const int atom_uc_cat = catom_array[atom].uc_category;
                        const int mat_uc_cat = create::internal::mp[mat].unit_cell_category;
								// check for within core shell range
								if(vmath::point_in_polygon_factor(x-x0,y-y0,factor, tmp_grain_pointx_array,tmp_grain_pointy_array,num_vertices)==true){
									if((cz>=minz) && (cz<maxz) && (atom_uc_cat == mat_uc_cat) ){
										catom_array[atom].include=true;
										catom_array[atom].material=mat;
										catom_array[atom].grain=grain;
									}
									// if set to clear atoms then remove atoms within radius
									else if(cs::fill_core_shell==false){
										catom_array[atom].include=false;
									}
								}
							}
						}
						// Check to see if site is within polygon
						else if(vmath::point_in_polygon_factor(x-x0,y-y0,1.0,tmp_grain_pointx_array,tmp_grain_pointy_array,num_vertices)==true){
							catom_array[atom].include=true;
							catom_array[atom].grain=grain;
						}
					}
				}
			}
		}
	}
	terminaltextcolor(GREEN);
	std::cout << "done!" << std::endl;
	terminaltextcolor(WHITE);
	zlog << "done!" << std::endl;

	// add final grain for continuous layer
	grain_coord_array.push_back(std::vector <double>());
	grain_coord_array[grain_coord_array.size()-1].push_back(0.0); // x
	grain_coord_array[grain_coord_array.size()-1].push_back(0.0); // y

	// check for continuous layer
	for(unsigned int atom=0; atom < catom_array.size(); atom++){
	  if(mp::material[catom_array[atom].material].continuous==true && catom_array[atom].include == false ){
	    catom_array[atom].include=true;
	    catom_array[atom].grain=int(grain_coord_array.size()-1);
	  }
	}

	// set number of grains
	grains::num_grains = int(grain_coord_array.size());

	// sort atoms by grain number
	create::internal::sort_atoms_by_grain(catom_array);

	return EXIT_SUCCESS;
}

} // End of cs namespace
