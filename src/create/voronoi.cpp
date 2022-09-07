//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2022. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <random>
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
#include "vmpi.hpp"
#include "vmath.hpp"
#include "voronoi.hpp"
#include <algorithm>

// create module headers
#include "internal.hpp"

namespace create_voronoi{
	bool parity=0;	/// left-right (0) or right-left (1) point initialisation
	bool rounded=false;
	double area_cutoff=0.8;
	double voronoi_sd=0.15;			/// Standard Deviation of voronoi grains
	bool include_boundary_grains_real = false;
}

namespace cs{

int voronoi_film(std::vector<cs::catom_t> & catom_array){

	// check calling of routine if error checking is activated
	if(err::check==true){
      terminaltextcolor(RED);
      std::cerr << "cs::voronoi_film has been called" << std::endl;
      terminaltextcolor(WHITE);
   }

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
   int grain = 0;
	// --------------------------------------------------------------------------
	// poisson distribution
	// --------------------------------------------------------------------------

	if (create::internal::grain_poission){

		std::ofstream file;
	   file.open("dist");

		int sdx = cs::system_dimensions[0];
		int sdy = cs::system_dimensions[1];

		//std::random_device rd;
      std::mt19937 gen(12345);
		//variance = exp ( 2.0 * mu + sigma * sigma ) * ( exp ( sigma * sigma ) - 1.0 );
		//mean = exp ( mu + 0.5 * sigma * sigma );
		std::lognormal_distribution<> d(log(grain_cell_size_x), grain_sd );
		std::cout << grain_cell_size_x << "\t" << grain_sd<< std::endl;
		double initial_grain_pos_x = sdx/2.0;
		double initial_grain_pos_y = sdy/2.0;
		double initial_grain_r = d(gen)/2.0;
		while (initial_grain_r > 2*grain_cell_size_x){
			initial_grain_r = d(gen)/2.0;
		}
		int maxattempts1 = 10;
		int maxattempts2 = 1000;

		std::vector <double> grains_x;
		std::vector <double> grains_y;
		std::vector <double> grains_r;
		std::vector <bool> active;
		grain=0;
		grains_x.push_back(initial_grain_pos_x);
		grains_y.push_back(initial_grain_pos_y);
		grains_r.push_back(initial_grain_r);
		active.push_back(true);
		grain_coord_array.push_back(std::vector <double>());
		grain_vertices_array.push_back(std::vector <std::vector <double> >());
		grain_coord_array[grain].push_back(initial_grain_pos_x);
		grain_coord_array[grain].push_back(initial_grain_pos_y);
		int num_active_grains = 1;
		grain++;
		double PI = 3.14159265;
		file << initial_grain_pos_x << '\t' << initial_grain_pos_y << '\t' << initial_grain_r << std::endl;
		for (int attempt = 0; attempt < maxattempts1; attempt ++){
			for (size_t g =0; g < grains_x.size() ; g++){
				double r = d(gen)/2.0;
			//	std::cout << r << "\t" << grain_cell_size_x+ grain_sd*grain_cell_size_x << "\t" << grain_cell_size_x << std::endl;
				while (r > grain_cell_size_x + grain_sd*grain_cell_size_x || r < grain_cell_size_x - grain_sd*grain_cell_size_x){
					r = d(gen)/2.0;
			//		std::cout << "A" <<std::endl;
 				}
			//	std::cout << r << std::endl;
				bool added =false;
				int attempt2 = 0;
				while (added == false && attempt2 < maxattempts2){
					attempt2 ++;
					int i = int(mtrandom::grnd()*grains_x.size());

				   double r1 = mtrandom::grnd();
			      double r2 = mtrandom::grnd();
			      double theta = r1*PI;
			      double phi = r2*PI*2;
					//file << number << std::endl;
					double min_distance = grains_r[i] + r;
					double max_distance = min_distance;
					double max_minus_min = max_distance - min_distance;
			      double rn = mtrandom::grnd();
			      double d = rn*max_minus_min;
			      double dx = sin(theta)*cos(phi);
			      double dy = sin(theta)*sin(phi);
			      double ddxdy = sqrt(dx*dx + dy*dy);
			      dx = dx/ddxdy;
			      dy = dy/ddxdy;
			      double x = grains_x[i] + d*dx + min_distance*dx;
			      double y = grains_y[i] + d*dy + min_distance*dy;
					if (x <= 2*sdx && y <= 2*sdy & x >= 0  - 2*grain_cell_size_x && y >= 0 - 2*grain_cell_size_y){
						int within =0;
				       for (size_t grain = 0; grain < grains_x.size(); grain ++ ){
				   		double dx2 = grains_x[grain] - x;
				   		double dy2 = grains_y[grain] - y;
				   		double dist = sqrt(dx2*dx2 + dy2*dy2);
				         if (dist < grains_r[grain] + r){
				   			within = 1;
				   			break;
				      	}
				   	}
						if (within ==0){
							bool xmove = true;
							bool ymove = true;
							double tempx = x;
							double tempy = y;
							while (xmove || ymove){
								if (x > sdx/2.0)	tempx = tempx - 1;
								else tempx = tempx + 1;
								for (size_t grain2 = 0; grain2 < grains_x.size(); grain2 ++ ){
								  double dx2 = grains_x[grain2] - tempx;
								  double dy2 = grains_y[grain2] - tempy;
								  double dist = sqrt(dx2*dx2 + dy2*dy2);
								  if (dist < grains_r[grain2] + r){
									  xmove = false;
									  break;
								  }
							  }
								if (y > sdy/2.0)	tempy = tempy - 1;
  								else tempy = tempy + 1;
  								for (size_t grain2 = 0; grain2 < grains_x.size(); grain2 ++ ){
  								  double dx2 = grains_x[grain2] - tempx;
  								  double dy2 = grains_y[grain2] - tempy;
  								  double dist = sqrt(dx2*dx2 + dy2*dy2);
  								  if (dist < grains_r[grain2] + r){
  									  ymove = false;
  									  break;
  								  }
							  }
							  if (xmove) x = tempx;
							  if (ymove) y = tempy;

							}

							//file << grain << '\t' << x/10 << '\t' << y/10 << '\t' << r/10 << std::endl;
					      grains_x.push_back(x);
							file << x << '\t' << y << '\t' << r << std::endl;
							grain_coord_array.push_back(std::vector <double>());
							grain_vertices_array.push_back(std::vector <std::vector <double> >());
							grain_coord_array[grain].push_back(x);
							grain_coord_array[grain].push_back(y);
					      grains_y.push_back(y);
					      grains_r.push_back(r);
					      active.push_back(true);
							num_active_grains ++;
							grain++;
							added = true;
						}
					}
				}
			}
		}
		double sumV = 0;
		double sumR = 0;
		for (size_t i = 0; i < grains_x.size(); i ++){
			double r = grains_r[i];
			double V = 3.14*r*r;
			sumV = sumV +V;
			sumR = sumR +r;
		}
		//------------------------------------------------------------------------
		// output voroni statistics to screen
		//------------------------------------------------------------------------
		double avR = sumR/grains_x.size();
		std::sort(grains_r.begin(),grains_r.end());
		int index = grains_r.size()/2;
		double Mr = (grains_r[index-1] + grains_r[index])/2;
		std::cout<< "Median grain radius: " << Mr << std::endl;
		std::cout<< "Mean grain radius: " << avR << std::endl;
		//------------------------------------------------------------------------
		double totalV = (2*sdx + 2*grain_cell_size_x)*(2*sdy + 2*grain_cell_size_y);
	 	double frac = sumV/totalV;
		std::cout<< " frac: " << frac << std::endl;
		std::cout << grain_cell_size_x/2.0 - avR << "\t" <<  2.0*(grain_cell_size_x/2.0 - avR)/grain_cell_size_x << std::endl;
		frac = frac + (grain_cell_size_x/2.0 - avR)/grain_cell_size_x;
	 	//	std::cout << sumV << '\t' << totalV << '\t' << frac  << "\t" << grains_x.size() << '\t'<<  grain_coord_array.size() << std::endl;

		for (size_t i = 0; i < grain_coord_array.size(); i ++){
		//	double r = grains_r[i];
		grain_coord_array[i][0] = grain_coord_array[i][0]* frac;
		grain_coord_array[i][1] = grain_coord_array[i][1]* frac;
		//	sd = sd + sqrt(r*r - Mr*Mr);
		//	file << grain << '\t' << grain_coord_array[i][0]/10 << '\t' << grain_coord_array[i][1]/10 << '\t' << grains_r[i]/10 << std::endl;
		}

	} // end of Poisson version

	else{
		//Calculate pointers
		for(int grain2=0;grain2<init_num_grains;grain2++){
			grain_coord_array.push_back(std::vector <double>());
			grain_coord_array[grain2].reserve(2);
			grain_vertices_array.push_back(std::vector <std::vector <double> >());
		}

		for (int x_particle=0;x_particle < num_x_particle;x_particle++){
			for (int y_particle=0;y_particle < num_y_particle;y_particle++){
				for (int particle_parity=0;particle_parity<2;particle_parity++){
					particle_coords[0] = (particle_parity)*delta_particle_x_parity + delta_particle_x*x_particle + vp*double(1-2*particle_parity)*delta_particle_x_parity;
					particle_coords[1] = (particle_parity)*delta_particle_y_parity + delta_particle_y*y_particle;
					grain_coord_array[grain].push_back(particle_coords[0]+grain_sd*mtrandom::gaussian()*delta_particle_x - cs::system_dimensions[0]);
					grain_coord_array[grain].push_back(particle_coords[1]+grain_sd*mtrandom::gaussian()*delta_particle_y - cs::system_dimensions[1]);
					grain++;
				}
			}
		}

		//-------------------------------------------------------
		// centre system on grain nearest the centre
		//-------------------------------------------------------
		// const double mpx = cs::system_dimensions[0] * 0.5;
		// const double mpy = cs::system_dimensions[1] * 0.5;
		// int nearest_grain = 0;
		// double nearest_distance = 1.0e99;
		// for(int grain = 0; grain < grain_coord_array[grain].size(); grain++){
		// 	double rx = grain_coord_array[grain][0];
		// 	double ry = grain_coord_array[grain][1];
		// 	double r2 = rx*rx + ry*ry;
		// 	// check for new nearest grain
		// 	if(r2 < nearest_distance){
		// 		nearest_distance = r2;
		// 		nearest_grain = grain;
		// 	}
		// }
		// // now shift all coordinates so that nearest grain is in the middle of the system
		// const double sx = mpx - grain_coord_array[nearest_grain][0];
		// const double sy = mpy - grain_coord_array[nearest_grain][1];
		// for(int grain = 0; grain < grain_coord_array[grain].size(); grain++){
		// 	grain_coord_array[grain][0] = grain_coord_array[grain][0]+sx;
		// 	grain_coord_array[grain][1] = grain_coord_array[grain][1]+sy;
		// }

	}



	// ----------------------- end of old version -------------------------------------

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
	create::internal::populate_vertex_points(grain_coord_array, grain_vertices_array, true);


	// Shrink Voronoi vertices in reduced coordinates to get spacing
	double shrink_factor = create::internal::voronoi_grain_size/(create::internal::voronoi_grain_size+create::internal::voronoi_grain_spacing);

	// Reduce vertices to relative coordinates
	for(unsigned int grain=0;grain<grain_coord_array.size();grain++){
		//grain_coord_array[grain][0]=0.0;
		//grain_coord_array[grain][1]=0.0;
	//	file1<< "h"<< grain << '\t' << grain_coord_array[grain][0] << '\t' << grain_coord_array[grain][1] << "\t" << delta_particle_x << std::endl;
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

	std::vector <double > R_med;

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

	std::cout <<"Generating Voronoi Grains" << std::flush;
	zlog << zTs() << "Generating Voronoi Grains" << std::flush;

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

// 	for(unsigned int grain=0;grain<grain_coord_array.size();grain++){
// 		file1<< "h"<< grain << '\t' << grain_coord_array[grain][0] << '\t' << grain_coord_array[grain][1] << "\t" << delta_particle_x << std::endl;
// }




	return EXIT_SUCCESS;
}

} // End of cs namespace
