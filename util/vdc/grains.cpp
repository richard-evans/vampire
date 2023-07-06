//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2023. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

// program header
#include "vdc.hpp"

//------------------------------------------------------------------------------
//	   Function to decide if test point is within polygon defined by points
//------------------------------------------------------------------------------
inline bool point_in_polygon(vdc::xy_t test, std::vector<vdc::xy_t>& points){  // list of points defining a polygon

	// Add tiny amount to atomic coordinates to include atoms at 0,0,0
	test.x += 1e-10;
	test.y += 1e-10;

	const int poly_sides = points.size();
	int j = poly_sides - 1 ;
	bool odd_nodes=false;

	for ( int i = 0 ; i < poly_sides; i++) {
		if ( ( points[i].y < test.y && points[j].y >= test.y ) || ( points[j].y < test.y && points[i].y >= test.y) ) {
			if (points[i].x+(test.y-points[i].y) / ( points[j].y - points[i].y ) * ( points[j].x - points[i].x) < test.x) {
				odd_nodes = !odd_nodes;
			}
		}
		j=i;
	}

  return odd_nodes;

}

namespace vdc{

   //---------------------------------------------------------------------------
   // Function to load grain shape file from VAMPIRE
   //---------------------------------------------------------------------------
   void load_grain_vertices(){

      //--------------------------------
      // set up input file stream
      //--------------------------------
      std::ifstream ifile;

      // open file and check for success
      ifile.open("grain_shapes.txt");
      // if failed inform user and exit
      if(!ifile.is_open()){
         std::cerr << "Error! grain_shapes.txt file not found required by vampire grain identification. Exiting" << std::endl;
         exit(234);
      }

      //--------------------------------
      // Load grain vertices
      //--------------------------------
      int grain_id = 0;
      std::string line;

      while(getline(ifile, line)){

         // create stringstream from line
         std::stringstream line_ss(line);
         int num_vertices = 0;
         line_ss >> num_vertices;

         // storage for vertices for this grain
         std::vector <xy_t> vertices;

         // loop over all vertices and read in xy pairs
         for(int v=0; v<num_vertices; v++){
            xy_t vv;
            line_ss >> vv.x >> vv.y;
            vertices.push_back(vv);
         }

         // save vertex list to grain vertices array
         vdc::grain_vertices_array.push_back(vertices);

      }

      // close input file
      ifile.close();

      return;

   }

   //---------------------------------------------------------------------------
   // Function to determine grain ID of each atom
   //---------------------------------------------------------------------------
   void determine_atom_grain_id(){

      // set last grain number (for atoms outside of grains)
      const int last_grain = vdc::grain_vertices_array.size()+1;

      // resize grain ID arrays
      vdc::grain.resize( vdc::type.size(), last_grain );
      vdc::nm_grain.resize( vdc::nm_type.size(), last_grain );

      for(size_t i=0; i < vdc::sliced_atoms_list.size(); i++){

         // get atom ID and coordinates
         const unsigned int atom = vdc::sliced_atoms_list[i];
         const double x = vdc::coordinates[3*atom + 0];
         const double y = vdc::coordinates[3*atom + 1];
         const double z = vdc::coordinates[3*atom + 2];
         const xy_t xy {x, y};

         // loop over all grains to see if atom is inside
         for(int g = 0; g < vdc::grain_vertices_array.size(); g++){
            if(point_in_polygon(xy, vdc::grain_vertices_array[g])){
               vdc::grain[atom] = g;
               break;
            }
         }

      }

      for(size_t i=0; i < vdc::sliced_nm_atoms_list.size(); i++){

         // get atom ID
         unsigned int atom = vdc::sliced_nm_atoms_list[i];

         const double x = vdc::nm_coordinates[3*atom + 0];
         const double y = vdc::nm_coordinates[3*atom + 1];
         const double z = vdc::nm_coordinates[3*atom + 2];
         const xy_t xy {x, y};

         // loop over all grains to see if atom is inside
         for(int g = 0; g < vdc::grain_vertices_array.size(); g++){
            if(point_in_polygon(xy, vdc::grain_vertices_array[g])){
               vdc::nm_grain[atom] = g;
               break;
            }
         }

      }

      return;
   }

   //---------------------------------------------------------------------------
   //---------------------------------------------------------------------------
   void generate_povray_grains(){

      return;
   }


}
