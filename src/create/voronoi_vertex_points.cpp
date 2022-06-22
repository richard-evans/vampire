//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>
#include <fstream>
#include <sstream>

// Vampire headers
#include "create.hpp"
#include "errors.hpp"
#include "qvoronoi.hpp"
#include "vio.hpp"
#include "voronoi.hpp"

// micromagnetic module headers
#include "internal.hpp"

namespace create{
namespace internal{

void populate_vertex_points(std::vector <std::vector <double> > & grain_coord_array,
                                std::vector <std::vector <std::vector <double> > > &  grain_vertices_array,
                                bool include_boundary_grains){ // controls is boundary grains are removed
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

	grain_file_sstr << "vor_grains.tmp";
	voronoi_file_sstr << "voronoi.tmp";

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

	//-----------------------------------------------------------------------------
	// Output grain coordindates and Scale to be unit length (0:1) on root process
	//-----------------------------------------------------------------------------

   bool root = false; // flag to indentify root process
   if(vmpi::my_rank == 0) root = true; // change flag to true on root process

   if(root){

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

   }

   //--------------------------------------------------------
   // Process Voronoi vertices from file
   //--------------------------------------------------------

	int dimensions,num_vertices,num_points,one;

   // stringstream stream declaration
   std::stringstream vertices_file;

   // fill input file stream with contents of file opened on master process
   vertices_file.str( vin::get_string(voronoi_filec, "voronoi", -1) );

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

   // scale and shift to centre of system
	for(int i=0;i<num_vertices;i++){
		vertex_array[i][0]=(vertex_array[i][0]+0.5)*scale_factor; // - cs::system_dimensions[0];
		vertex_array[i][1]=(vertex_array[i][1]+0.5)*scale_factor; // - cs::system_dimensions[0];
	}

   //--------------------------------------
   // Read in Voronoi vertex associations
   //--------------------------------------
   for(int i=0;i<num_grains;i++){
      int num_assoc_vertices;	// Number of vertices associated with point i
      int vertex_number;		// temporary vertex number
      vertices_file >> num_assoc_vertices;
      bool inf=false;

      //std::cout << i << '\t' << num_assoc_vertices <<std::endl;
      for(int j=0;j<num_assoc_vertices;j++){
         vertices_file >> vertex_number;
         //grain_vertices_array[i].push_back(std::vector <double>());
         //grain_vertices_array[i][j].push_back(vertex_array[vertex_number][0]);
         //grain_vertices_array[i][j].push_back(vertex_array[vertex_number][1]);
         // check for unbounded grains
         if(vertex_number==0) inf=true;
         //std::cout <<i << '\t' <<  vertex_array[vertex_number][0] << '\t' << vertex_array[vertex_number][1] << '\t' << cs::system_dimensions[0] << '\t' << cs::system_dimensions[1] <<std::endl;
         // Check for whether boundary grains are included or removed
         if(create_voronoi::include_boundary_grains_real){
            //std::cout << "A" << std::endl;
            // check for bounded grains with vertices outside bounding box and truncate to system edges
            if(vertex_array[vertex_number][0] < -0.5*cs::system_dimensions[0]) vertex_array[vertex_number][0] = 0.0;
            if(vertex_array[vertex_number][0] >  1.5*cs::system_dimensions[0]) vertex_array[vertex_number][0] = cs::system_dimensions[0];
            if(vertex_array[vertex_number][1] < -0.5*cs::system_dimensions[1]) vertex_array[vertex_number][1] = 0.0;
            if(vertex_array[vertex_number][1] >  1.5*cs::system_dimensions[1]) vertex_array[vertex_number][1] = cs::system_dimensions[1];
            // check for bounded grains with vertices outside bounding box
            if(vertex_array[vertex_number][0]<0.0) {
               vertex_array[vertex_number][0] = 0.0;
            }
            if (vertex_array[vertex_number][0]>cs::system_dimensions[0]){
               vertex_array[vertex_number][0] = cs::system_dimensions[0];
            }
            if(vertex_array[vertex_number][1]<0.0){
               vertex_array[vertex_number][1] = 0.0;
            }
            if (vertex_array[vertex_number][1]>cs::system_dimensions[1]) {
               vertex_array[vertex_number][1] = cs::system_dimensions[1];
            }

            // but only include grains with a centre of mass inside the system dimensions
            const double tol = create::internal::voronoi_grain_size * 0.25; // allow grains with half the mean grain size
            bool inx = grain_coord_array[i][0] >= -tol && grain_coord_array[i][0] < cs::system_dimensions[0] + tol;
            bool iny = grain_coord_array[i][1] >= -tol && grain_coord_array[i][1] < cs::system_dimensions[1] + tol;
            if( !(inx && iny) ) inf = true;

         }
          else{
             // check for bounded grains with vertices outside bounding box
             if((vertex_array[vertex_number][0]<0.0) || (vertex_array[vertex_number][0]>cs::system_dimensions[0])) inf=true;
             if((vertex_array[vertex_number][1]<0.0) || (vertex_array[vertex_number][1]>cs::system_dimensions[1])) inf=true;
          }
         grain_vertices_array[i].push_back(std::vector <double>());
         grain_vertices_array[i][j].push_back(vertex_array[vertex_number][0]);
         grain_vertices_array[i][j].push_back(vertex_array[vertex_number][1]);
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
   if(root){
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
		  std::cerr << "Error removing temporary grain file" << std::endl;
		  terminaltextcolor(WHITE);
	   }
   }
   if(root){
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
		  std::cerr << "Error removing temporary voronoi file" << std::endl;
		  terminaltextcolor(WHITE);
	   }
   }

   // Recalculate grain coordinates as average of vertices
   for(unsigned int grain=0;grain<grain_coord_array.size();grain++){
      grain_coord_array[grain][0]=0.0; // could be a temporary to avoid multiple array writes?
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

   return;
}

} // end of namespace internal
} // end of namespace create
