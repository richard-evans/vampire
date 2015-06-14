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
#include "grains.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "random.hpp"
#include "vmpi.hpp"
#include "vmath.hpp"
#include "vio.hpp"
#include "qvoronoi.hpp"


#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>



namespace create_voronoi{
	bool parity=0;	/// left-right (0) or right-left (1) point initialisation
	bool rounded=false;
	double area_cutoff=0.8;
	double voronoi_sd=0.15;			/// Standard Deviation of voronoi grains
	
}

namespace cs{
	
   //----------------------------------------
   // Function prototypes
   //----------------------------------------
   int populate_vertex_points(std::vector <std::vector <double> > &, std::vector <std::vector <std::vector <double> > > &);

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
	double size = cs::particle_scale+cs::particle_spacing;
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

	//-----------------------------------------------
	// Calculate Voronoi construction using qhull
	//-----------------------------------------------
	//populate_vertex_points(
	//	grains::num_grains,init_grain_coord_array,init_grain_pointx_array,
	//	init_grain_pointy_array,init_num_assoc_vertices_array,vbox_dimensions);
	
	populate_vertex_points(grain_coord_array, grain_vertices_array);

	// Recalculate grain coordinates as average of vertices
	for(unsigned int grain=0;grain<grain_coord_array.size();grain++){
		grain_coord_array[grain][0]=0.0;
		grain_coord_array[grain][1]=0.0;
		const int nv = grain_vertices_array[grain].size();
		// Exclude grains with zero vertices
		if(nv!=0){
			//std::cout << grain_vertices_array[grain].size() << " " << grain_vertices_array[grain].capacity() << std::endl;
			for(int vertex=0;vertex<nv;vertex++){
			//	std::cout << grain_vertices_array[grain][vertex][0] << "\t" << grain_vertices_array[grain][vertex][1] << std::endl;
				grain_coord_array[grain][0]+=grain_vertices_array[grain][vertex][0];
				grain_coord_array[grain][1]+=grain_vertices_array[grain][vertex][1];
			}
			grain_coord_array[grain][0]/=double(nv);
			grain_coord_array[grain][1]/=double(nv);
		}
	}

	// Shrink Voronoi vertices in reduced coordinates to get spacing
	double shrink_factor = cs::particle_scale/(cs::particle_scale+cs::particle_spacing);
	
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

	double tmp_grain_pointx_array[max_vertices];
	double tmp_grain_pointy_array[max_vertices];

	// calculate grain rounding
	if(create_voronoi::rounded==true){
		for(unsigned int grain=0;grain<grain_coord_array.size();grain++){

			// get number of vertices for each grain
			const int nv = grain_vertices_array[grain].size();

			// Exclude grains with zero vertices
			if(nv!=0){

				// allocate temporary array for area calculation 
				std::vector<std::vector <double> > rnd;
				rnd.resize(48);
				for(int idx=0; idx<48;idx++){
					rnd[idx].resize(2,0.0);
				}
				
				double area_frac=create_voronoi::area_cutoff;
				double deltar=0.5; // Angstroms
				double radius=0.0; // Strating radius
				double area=0.0; // Starting area

				// Set temporary vertex coordinates
				int num_vertices = grain_vertices_array[grain].size();
					for(int vertex=0;vertex<num_vertices;vertex++){
						tmp_grain_pointx_array[vertex]=grain_vertices_array[grain][vertex][0]; //-grain_coord_array[grain][0];
						tmp_grain_pointy_array[vertex]=grain_vertices_array[grain][vertex][1]; //-grain_coord_array[grain][1];
					}

				// calculate voronoi area
				double varea=0.0;
				for(int r=0;r<1000;r++){
					radius+=deltar;

					for(int i=0;i<48;i++){
						double theta = 2.0*M_PI*double(i)/48.0;
						double x = radius*cos(theta);
						double y = radius*sin(theta);

						// Check to see if site is within polygon
						if(vmath::point_in_polygon(x,y,tmp_grain_pointx_array,tmp_grain_pointy_array,num_vertices)==true){
							rnd[i][0]=x;
							rnd[i][1]=y;
						}
					}
					
					//update area
					varea=0.0;
					for(int i=0;i<48;i++){
						int nvi = i+1;
						if(nvi>=48) nvi=0;
						varea+=0.5*sqrt((rnd[nvi][0]-rnd[i][0])*(rnd[nvi][0]-rnd[i][0]))*sqrt((rnd[nvi][1]-rnd[i][1])*(rnd[nvi][1]-rnd[i][1]));
					}
				}

				// reset polygon positions and radius
				for(int idx=0; idx<48;idx++){
					rnd[idx][0]=0.0;
					rnd[idx][1]=0.0;
				}
				radius=0.0;
				// expand polygon
				for(int r=0;r<100;r++){
				if(area<area_frac*varea){
					radius+=deltar;
					
					//loop over coordinates
					for(int i=0;i<48;i++){
						double theta = 2.0*M_PI*double(i)/48.0;
						double x = radius*cos(theta);
						double y = radius*sin(theta);

						// Check to see if site is within polygon
						if(vmath::point_in_polygon(x,y,tmp_grain_pointx_array,tmp_grain_pointy_array,num_vertices)==true){
							rnd[i][0]=x;
							rnd[i][1]=y;
						}
					}
					
					//update area
					area=0.0;
					for(int i=0;i<48;i++){
						int nvi = i+1;
						if(nvi>=48) nvi=0;
						area+=0.5*sqrt((rnd[nvi][0]-rnd[i][0])*(rnd[nvi][0]-rnd[i][0]))*sqrt((rnd[nvi][1]-rnd[i][1])*(rnd[nvi][1]-rnd[i][1]));
					}
				}
			}

			// set new polygon points
			grain_vertices_array[grain]=rnd;
				
			} // end of check for nv=0
		
		} // end of grain loop
	} // end of rounding if

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
	
	std::cout <<"Generating Voronoi Grains";
	zlog << zTs() << "Generating Voronoi Grains";

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
				tmp_grain_pointx_array[vertex]=grain_vertices_array[grain][vertex][0]+grain_coord_array[grain][0];
				tmp_grain_pointy_array[vertex]=grain_vertices_array[grain][vertex][1]+grain_coord_array[grain][1];
				int x = int(tmp_grain_pointx_array[vertex]/unit_cell.dimensions[0]);
				int y = int(tmp_grain_pointy_array[vertex]/unit_cell.dimensions[1]);
				if(x < minx) minx = x;
				if(x > maxx) maxx = x;
				if(y < miny) miny = y;
				if(y > maxy) maxy = y;
			}

			// loopover cells
			for(int i=minx;i<=maxx;i++){
				for(int j=miny;j<=maxy;j++){
				
					// loop over atoms in cells;
					for(int id=0;id<supercell_array[i][j].size();id++){
						int atom = supercell_array[i][j][id];
						
						double x = catom_array[atom].x;
						double y = catom_array[atom].y;
						
						// Check to see if site is within polygon
						if(vmath::point_in_polygon(x,y,tmp_grain_pointx_array,tmp_grain_pointy_array,num_vertices)==true){
							catom_array[atom].include=true;
							catom_array[atom].grain=grain;
						}
						// Check for continuous media
						/*else if(mp::material[catom_array[atom].material].continuous==true){
							const int geo=mp::material[catom_array[atom].material].geometry;

							if(geo==0){
								catom_array[atom].include=true;
							}
							else{
								double x = catom_array[atom].x;
								double y = catom_array[atom].y;
								double px[50];
								double py[50];
								// Initialise polygon points
								for(int p=0;p<geo;p++){
									px[p]=mp::material[catom_array[atom].material].geometry_coords[p][0]*cs::system_dimensions[0];
									py[p]=mp::material[catom_array[atom].material].geometry_coords[p][1]*cs::system_dimensions[1];
								}
								if(vmath::point_in_polygon(x,y,px,py,geo)==true){
									catom_array[atom].include=true;
									catom_array[atom].grain=grain;
								}
							}
							//catom_array[atom].include=true;
							catom_array[atom].grain=grain;
							}*/
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
	for(int atom=0; atom < catom_array.size(); atom++){
	  if(mp::material[catom_array[atom].material].continuous==true && catom_array[atom].include == false ){
	    catom_array[atom].include=true;
	    catom_array[atom].grain=int(grain_coord_array.size()-1);
	  }
	}
  
	// set number of grains
	grains::num_grains = int(grain_coord_array.size());

	// sort atoms by grain number
	sort_atoms_by_grain(catom_array);

	return EXIT_SUCCESS;	
}

int populate_vertex_points(std::vector <std::vector <double> > & grain_coord_array, std::vector <std::vector <std::vector <double> > > &  grain_vertices_array){
	//========================================================================================================
	//		 				Function to populate voronoi vertices for grains using qhull
	//
	//														Version 1.0
	//
	//												R F Evans 15/07/2009
	//
	//========================================================================================================
	//		Locally allocated variables: vertex_array
	//========================================================================================================

	const int num_grains=grain_coord_array.size();
	std::stringstream grain_file_sstr;
	std::stringstream voronoi_file_sstr;

	grain_file_sstr << "grains." << vmpi::my_rank << ".tmp";
	voronoi_file_sstr << "voronoi." << vmpi::my_rank << ".tmp";

	string grain_file = grain_file_sstr.str();
	string voronoi_file = voronoi_file_sstr.str();
	const char* grain_filec = grain_file.c_str();
	const char* voronoi_filec = voronoi_file.c_str();

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){
		terminaltextcolor(RED);
		std::cerr << "cs::populate_vertex_points has been called" << std::endl;
		terminaltextcolor(WHITE);
	}

	//--------------------------------------
	// Scale grain coordinates
	//--------------------------------------

	double scale_factor=0.0;

	//--------------------------------------------------------------
	// Calculate max grain coordinate
	//--------------------------------------------------------------
	for(int i=0;i<num_grains;i++){
		if(grain_coord_array[i][0]>scale_factor){
			scale_factor=grain_coord_array[i][0];
		}
		if(grain_coord_array[i][1]>scale_factor){
			scale_factor=grain_coord_array[i][1];
		}
	}

	//--------------------------------------------------------------
	// Output grain coordindates and Scale to be unit length (0:1)
	//--------------------------------------------------------------

	std::ofstream grain_coord_file;
  	grain_coord_file.open (grain_filec);

	grain_coord_file << "#============================================" << std::endl;
	grain_coord_file << "# Grain Coordinate file for input into qhull" << std::endl;
	grain_coord_file << "#============================================" << std::endl;
	grain_coord_file << "# Grain Scaling Factor" << std::endl;
	grain_coord_file << "#\t" << scale_factor << std::endl;
	grain_coord_file << 2 << std::endl;
	grain_coord_file << "# Number of Grains" << std::endl;
	grain_coord_file << num_grains << std::endl;
	grain_coord_file << "# Grains Coordinates" << std::endl;

	for(int i=0;i<num_grains;i++){
		grain_coord_file << grain_coord_array[i][0]/scale_factor-0.5 << "\t";
 		grain_coord_file << grain_coord_array[i][1]/scale_factor-0.5 << std::endl;
	}

	grain_coord_file.close();

   //--------------------------------------------------------
   // Use voronoi library creating an import and export temporary files
   //--------------------------------------------------------
   FILE *inputqv, *outputqv;
   int qargc=3;
   const char *qargv[3]={"qvoronoi", "-o", "-Fv"};
   inputqv=fopen(grain_file.c_str(),"r");
   outputqv=fopen(voronoi_file.c_str(),"w");
   qvoronoi(qargc, const_cast<char**>(qargv), inputqv, outputqv);
   fclose(outputqv);
   fclose(inputqv);

   //--------------------------------------------------------
   // Read in number of Voronoi vertices
   //--------------------------------------------------------

	int dimensions,num_vertices,num_points,one;

	std::ifstream vertices_file;
  	vertices_file.open (voronoi_filec);
	vertices_file >> dimensions;
	vertices_file >> num_vertices >> num_points >> one;


	//----------------------------------------------------------
	// Allocate vertex_array
	//----------------------------------------------------------
	std::vector<std::vector<double> > vertex_array;

    	vertex_array.resize(num_vertices);
	for(int i=0; i<num_vertices; i++) vertex_array[i].resize(2);

	//--------------------------------------
	// Read in Voronoi vertices and rescale
	//--------------------------------------

	for(int i=0;i<num_vertices;i++){
		vertices_file >> vertex_array[i][0];
		vertices_file >> vertex_array[i][1];
	}

	for(int i=0;i<num_vertices;i++){
		vertex_array[i][0]=(vertex_array[i][0]+0.5)*scale_factor;
		vertex_array[i][1]=(vertex_array[i][1]+0.5)*scale_factor;
	}

   //--------------------------------------
   // Read in Voronoi vertex associations
   //--------------------------------------
   for(int i=0;i<num_grains;i++){
      int num_assoc_vertices;	// Number of vertices associated with point i
      int vertex_number;		// temporary vertex number
      vertices_file >> num_assoc_vertices;
      bool inf=false;
      for(int j=0;j<num_assoc_vertices;j++){
         vertices_file >> vertex_number;
         grain_vertices_array[i].push_back(std::vector <double>());
         grain_vertices_array[i][j].push_back(vertex_array[vertex_number][0]);
         grain_vertices_array[i][j].push_back(vertex_array[vertex_number][1]);
         // check for unbounded grains
         if(vertex_number==0) inf=true;
         // check for bounded grains with vertices outside bounding box
         if((vertex_array[vertex_number][0]<0.0) || (vertex_array[vertex_number][0]>cs::system_dimensions[0])) inf=true;
         if((vertex_array[vertex_number][1]<0.0) || (vertex_array[vertex_number][1]>cs::system_dimensions[1])) inf=true;
      }

      //-------------------------------------------------------------------
      // Set unbounded grains to zero vertices for later removal
      //-------------------------------------------------------------------
      if(inf==true){
         for(int k=0;k<num_assoc_vertices;k++){
            grain_vertices_array[i][k].resize(0);
         }
         grain_vertices_array[i].resize(0);
      }
   }
   vertices_file.close();
   //-------------------------------------------------------------------
   // Remove Temporary voronoi files
   //-------------------------------------------------------------------
   {
      std::stringstream rmfsstr;
      #ifdef WIN_COMPILE
         //rmfsstr << "del /f " << grain_file;
      #else
         rmfsstr << "rm -f " << grain_file;
      #endif
      string rmfstr = rmfsstr.str();
      const char* rmfcstr = rmfstr.c_str();
      int sysstat=system(rmfcstr);
      if(sysstat!=0) {
		  terminaltextcolor(RED);
		  std::cerr << "Error removing temporary grain files" << std::endl;
		  terminaltextcolor(WHITE);
	  }
   }
   {
      std::stringstream rmfsstr;
      #ifdef WIN_COMPILE
         rmfsstr << "del /f " << voronoi_file;
      #else
         rmfsstr << "rm -f " << voronoi_file;
      #endif
      string rmfstr = rmfsstr.str();
      const char* rmfcstr = rmfstr.c_str();
      int sysstat=system(rmfcstr);
      if(sysstat!=0) {
		  terminaltextcolor(RED);
		  std::cerr << "Error removing temporary voronoi files" << std::endl;
		  terminaltextcolor(WHITE);
	  }
   }
   return EXIT_SUCCESS;

}



} // End of cs namespace

