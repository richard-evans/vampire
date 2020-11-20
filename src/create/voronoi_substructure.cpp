//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2016. All rights reserved.
//
//-----------------------------------------------------------------------------

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
#include "qvoronoi.hpp"
#include "random.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"
#include "voronoi.hpp"

// Internal create header
#include "internal.hpp"

namespace create{
namespace internal{

bool in_pill(double x, double y, double z, double px, double py, double ptz, double pbz, double r);

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

	//---------------------------------------------------
	// Local constants
	//---------------------------------------------------
	const int max_vertices=50;
	double grain_sd=create_voronoi::voronoi_sd;

	// Set number of particles in x and y directions
	double size = create::internal::voronoi_grain_size + create::internal::voronoi_grain_spacing;
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

	// Calculate Voronoi construction using qhull not including boundary grains (false) - to be fixed!
	create::internal::populate_vertex_points(grain_coord_array, grain_vertices_array, false);

	// Shrink Voronoi vertices in reduced coordinates to get spacing
	double shrink_factor = create::internal::voronoi_grain_size/(create::internal::voronoi_grain_size+create::internal::voronoi_grain_spacing);

	// Reduce vertices to relative coordinates
	for(unsigned int grain=0;grain<grain_coord_array.size();grain++){
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
		int cx = int (catom_array[atom].x/cs::unit_cell.dimensions[0]);
		int cy = int (catom_array[atom].y/cs::unit_cell.dimensions[1]);
		supercell_array.at(cx).at(cy).push_back(atom);
	}

	// Determine order for core-shell grains
   std::list<core_radius_t> material_order(0);
   for(int mat=0;mat<mp::num_materials;mat++){
      core_radius_t tmp;
      tmp.mat=mat;
      tmp.radius=mp::material[mat].core_shell_size;
      material_order.push_back(tmp);
   }
   // sort by increasing radius
   material_order.sort(compare_radius);

	std::cout <<"Generating voronoi substructure";
	zlog << zTs() << "Generating voronoi substructure";

   // arrays to store list of grain vertices
   double tmp_grain_pointx_array[max_vertices];
	double tmp_grain_pointy_array[max_vertices];

   // array to store if atoms are included in substructure (assume not)
   std::vector<bool> insub(catom_array.size(),false);

   const double ssz = cs::system_dimensions[2];

   //------------------------------------------------------------------------------
   // Set 3D structure for grains
   //------------------------------------------------------------------------------
   const double sphere_radius = create::internal::voronoi_grain_substructure_crystallization_radius; //*create::internal::voronoi_grain_size*radius_factor;

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
				int x = int((tmp_grain_pointx_array[vertex]+grain_coord_array[grain][0])/cs::unit_cell.dimensions[0]);
				int y = int((tmp_grain_pointy_array[vertex]+grain_coord_array[grain][1])/cs::unit_cell.dimensions[1]);
				if(x < minx) minx = x;
				if(x > maxx) maxx = x;
				if(y < miny) miny = y;
				if(y > maxy) maxy = y;
			}

			// determine coordinate offset for grains
			const double x0 = grain_coord_array[grain][0];
			const double y0 = grain_coord_array[grain][1];

         // const double sphere_x = grain_coord_array[grain][0]; // unused variable
         // const double sphere_y = grain_coord_array[grain][1]; // unused variable

			// copy overlap of substructure grains to a local constant
			const double overlap = create::internal::voronoi_grain_substructure_overlap_factor;

			// loopover cells
			for(int i=minx;i<=maxx;i++){
				for(int j=miny;j<=maxy;j++){

					// loop over atoms in cells and z
               for(unsigned int id=0;id<supercell_array[i][j].size();id++){
                  const int atom = supercell_array[i][j][id];

                  // Get atomic position
                  const double x = catom_array[atom].x;
                  const double y = catom_array[atom].y;
                  const double z = catom_array[atom].z;
                  const double frh = z/ssz;
                  const double nucleation_height = create::internal::mp[catom_array[atom].material].voronoi_grain_substructure_nucleation_height;

                  // Remove core-shell code in structure for now. This means that core shell
                  // structures can still be applied to the superstructure, eg a dot or particle
                  /*if(mp::material[catom_array[atom].material].core_shell_size>0.0){
                     // Iterate over materials
                     for(std::list<core_radius_t>::iterator it = material_order.begin(); it !=  material_order.end(); it++){
                        int mat = (it)->mat;
                        double factor = mp::material[mat].core_shell_size;
                        double maxz=create::internal::mp[mat].max*cs::system_dimensions[2];
                        double minz=create::internal::mp[mat].min*cs::system_dimensions[2];
                        double cz=catom_array[atom].z;
                        // calculate reduced ranges for materials with small offset to prevent dangling atoms
                        double rminz = create::internal::mp[mat].min-0.01;
                        double rmaxz = create::internal::mp[mat].max+0.01;
                        double factor_radius = 0.0;
                        if(frh > nucleation_height){
                           // multiply by small factor to ensure grains touch at boundary for zero spacing
                           factor_radius = 1.04*pow(1.0+((nucleation_height-frh)/(rmaxz-nucleation_height)),sphere_radius);
                        }
                        else{
                           factor_radius = 1.04*pow((1.0-(frh-nucleation_height)/(rminz-nucleation_height)),sphere_radius);
                        }
                        // check for within core shell range
                        if(vmath::point_in_polygon_factor(x-x0,y-y0,factor*factor_radius, tmp_grain_pointx_array,tmp_grain_pointy_array,num_vertices)){ //}; // && in_pill(x, y, z, sphere_x, sphere_y, top_sphere_z, bot_sphere_z, factor*sphere_radius)){
                        if((cz>=minz) && (cz<maxz)){
                           insub[atom] = true;
                           catom_array[atom].material=mat;
                        }
                        // if set to clear atoms then remove atoms within radius
                        else if(cs::fill_core_shell==false){
                           insub[atom] = false;
                        }
                     }
                  }
                  }*/
                  // Check to see if site is within polygon
                  //else{
                  int mat = catom_array[atom].material;
                  //double maxz=create::internal::mp[mat].max*cs::system_dimensions[2]; // unused variable
                  //double minz=create::internal::mp[mat].min*cs::system_dimensions[2]; // unused variable
                  // calculate reduced ranges for materials with small offset to prevent dangling atoms
                  double rminz = create::internal::mp[mat].min-0.01;
                  double rmaxz = create::internal::mp[mat].max+0.01;
                  double factor_radius = 0.0;
                  if(frh > nucleation_height){
                     // multiply by small factor to ensure grains touch at boundary for zero spacing
                     factor_radius = 1.04*pow(1.0+((nucleation_height-frh)/(rmaxz-nucleation_height)),sphere_radius);
                  }
                  else{
                     factor_radius = 1.04*pow((1.0-(frh-nucleation_height)/(rminz-nucleation_height)),sphere_radius);
                  }
                  if(vmath::point_in_polygon_factor(x-x0,y-y0,overlap*factor_radius,tmp_grain_pointx_array,tmp_grain_pointy_array,num_vertices)){
                     insub[atom] = true;
                  }
                  //}
					}
				}
			}
		}
	}

	terminaltextcolor(GREEN);
	std::cout << "done!" << std::endl;
	terminaltextcolor(WHITE);
	zlog << "done!" << std::endl;

   // Now fill in with fill materials
   for(int mat=0;mat<mp::num_materials;mat++){
      if(create::internal::mp[mat].sub_fill){
         double min = create::internal::mp[mat].min*cs::system_dimensions[2];
         double max = create::internal::mp[mat].max*cs::system_dimensions[2];

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

bool in_pill(double x, double y, double z, double px, double py, double ptz, double pbz, double r){

   const double r2 = r*r;
   const double dx2 = (x - px)*(x - px);
   const double dy2 = (y - py)*(y - py);
   const double dtz2 = (z - ptz)*(z - ptz);
   const double dbz2 = (z - pbz)*(z - pbz);

   //bool in_r = dx2 + dy2 < r2; // unused variable

   if(r2 == false) return false;
   bool in_cylinder = z < ptz && z > pbz;
   bool in_top = dx2 + dy2 + dtz2 < r2;
   bool in_bot = dx2 + dy2 + dbz2 < r2;
   if(in_cylinder || in_top || in_bot) return true;
   return false;
}

} // end of namespace internal
} // end of namespace create
