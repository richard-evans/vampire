#include "create.hpp"
#include "grains.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "random.hpp"
#include "vmpi.hpp"
#include "vmath.hpp"

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

namespace create_voronoi{
	bool parity=0;	// left-right (0) or right-left (1) point initialisation
	double voronoi_sd=0.15;			// Standard Deviation of voronoi grains
}

namespace cs{
	
int voronoi_film(std::vector<cs::catom_t> & catom_array){
	
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "cs::voronoi_film has been called" << std::endl;}	
	
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

	//----------------------------------------
	// Function prototypes
	//----------------------------------------

	int populate_vertex_points(std::vector <std::vector <double> > &, std::vector <std::vector <std::vector <double> > > &);

	//---------------------------------------------------
	// Local constants
	//---------------------------------------------------
	const int max_vertices=50;
	double grain_sd=create_voronoi::voronoi_sd;

	// Set number of particles in x and y directions
	double size = material_parameters::particle_scale+material_parameters::particle_spacing;
	double grain_cell_size_x = size;
	double grain_cell_size_y = sqrt(3.0)*size;

	int num_x_particle = 4+iround(mp::system_dimensions[0]/(grain_cell_size_x));
	int num_y_particle = 4+iround(mp::system_dimensions[1]/(grain_cell_size_y));

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
		std::cerr << "Error! - No grains found in structure - Increase system dimensions" << std::endl;
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

	//std::cout << "here" << std::endl;
	
	// Shrink Voronoi vertices in reduced coordinates to get spacing
	double shrink_factor = mp::particle_scale/(mp::particle_scale+mp::particle_spacing);
	
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

	std::cout <<"Generating Voronoi Grains";

	double tmp_grain_pointx_array[max_vertices];
	double tmp_grain_pointy_array[max_vertices];
	const int num_atoms = catom_array.size();	
	int nearest_grain=0;
	
	for(unsigned int atom=0;atom<catom_array.size();atom++){
			if((atom%(num_atoms/10))==0){
				std::cout << "." << std::flush;
			}

		// Loop over grains to find if point is in grain
		for(unsigned int grain=0;grain<grain_coord_array.size();grain++){
			// Exclude grains with zero vertices
			if(grain_vertices_array[grain].size()!=0){
				double x = catom_array[atom].x-grain_coord_array[grain][0];
				double y = catom_array[atom].y-grain_coord_array[grain][1];

				// Set temporary vertex coordinates
				int num_vertices = grain_vertices_array[grain].size();
				for(int vertex=0;vertex<num_vertices;vertex++){
					tmp_grain_pointx_array[vertex]=grain_vertices_array[grain][vertex][0];
					tmp_grain_pointy_array[vertex]=grain_vertices_array[grain][vertex][1];
				}

				// Check to see if site is within polygon
				if(vmath::point_in_polygon(x,y,tmp_grain_pointx_array,tmp_grain_pointy_array,num_vertices)==true){
					catom_array[atom].include=true;
					catom_array[atom].grain=grain;
					break;
				}
				// Check for continuous media
				else if(mp::material[catom_array[atom].material].continuous==true){
					double new_range = x*x+y*y;
					double ox = catom_array[atom].x-grain_coord_array[nearest_grain][0];
					double oy = catom_array[atom].y-grain_coord_array[nearest_grain][1];
					double old_range = ox*ox+oy*oy;
					if(new_range<=old_range){
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
					      px[p]=mp::material[catom_array[atom].material].geometry_coords[p][0]*mp::system_dimensions[0];
					      py[p]=mp::material[catom_array[atom].material].geometry_coords[p][1]*mp::system_dimensions[1];
					    }
					    if(vmath::point_in_polygon(x,y,px,py,geo)==true){
					      catom_array[atom].include=true;
					      catom_array[atom].grain=grain;
					    }
					  }
						//catom_array[atom].include=true;
						catom_array[atom].grain=grain;
						nearest_grain=grain;
					}
				}
			}
		}
		// Identify atoms by grain parity
		//if(catom_array[atom].material==2){
		//	catom_array[atom].material=catom_array[atom].grain%2+2;
		//}
	}

	std::cout << "done!" << std::endl;

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
	if(err::check==true){std::cout << "cs::populate_vertex_points has been called" << std::endl;}	

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

	//--------------------------------------
	// Test for qhull installation
	// (Except for parallel version)
	//--------------------------------------
	#ifdef MPICF
	#else
	int test_qhull;

	test_qhull=system(VORONOI" 1> /dev/null 2> /dev/null");
	if(test_qhull!=0){
		std::cerr << "Error! - qvoronoi does not seem to be installed, exiting" << std::endl; 
		err::vexit();
	}
	#endif

	//--------------------------------------
	// Generate voronoi vertices with qhull
	//--------------------------------------

	std::stringstream vsstr;
	vsstr << "cat " << grain_file << " | "VORONOI" -o -Fv > " << voronoi_file;
	string vstr = vsstr.str();
	const char* vcstr = vstr.c_str();
	//std::cout << vcstr << std::endl;

	//system("cat grains.tmp | qvoronoi -o -Fv > voronoi.tmp");
	int sysstat = system(vcstr);
	if(sysstat!=0){
	  std::cerr << "Error initiating qvoronoi, exiting" << std::endl;
	  err::vexit();
	}
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
	double** vertex_array;

	try{vertex_array=new double*[num_vertices];
    	for(int i=0; i<num_vertices; i++)vertex_array[i]=new double[2];}
  	catch(...){std::cerr << "error allocating vertex_array" << std::endl;err::vexit();}

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
			if((vertex_array[vertex_number][0]<0.0) || (vertex_array[vertex_number][0]>mp::system_dimensions[0])) inf=true;
			if((vertex_array[vertex_number][1]<0.0) || (vertex_array[vertex_number][1]>mp::system_dimensions[1])) inf=true;
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

	//-------------------------------------------------------------------
	// Remove Temporary voronoi files
	//-------------------------------------------------------------------
	{
	std::stringstream rmfsstr;
	rmfsstr << "rm -f " << grain_file;
	string rmfstr = rmfsstr.str();
	const char* rmfcstr = rmfstr.c_str();
	int sysstat=system(rmfcstr);
	if(sysstat!=0) std::cerr << "Error removing temporary grain files" << std::endl;
	}
	{
	std::stringstream rmfsstr;
	rmfsstr << "rm -f " << voronoi_file;
	string rmfstr = rmfsstr.str();
	const char* rmfcstr = rmfstr.c_str();
	int sysstat=system(rmfcstr);
	if(sysstat!=0) std::cerr << "Error removing temporary voronoi files" << std::endl;
	}
		
	//system("rm -f grains.tmp");
	//system("rm -f voronoi.tmp");

	//-------------------------------------------------------------------
	// Deallocate vertex array
	//-------------------------------------------------------------------
	
 	try{for(int i=0;i<num_vertices;i++) delete [] vertex_array[i];
    	delete [] vertex_array;
    	vertex_array=NULL;
   }
  	catch(...){std::cerr << "error deallocating vertex_array" << std::endl;err::vexit();}
	
	return EXIT_SUCCESS;

}



} // End of cs namespace

