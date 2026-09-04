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

// Vampire headers
#include "grains.hpp"
#include "vmpi.hpp"
#include "vio.hpp"

// grains module headers
#include "internal.hpp"

//------------------------------------------------------------------------------
// Function to write grain vertices to a file in the format
// num_vertices vx1 vy1 vx2 vy2 ...
//------------------------------------------------------------------------------
// Expects the vertices array in the format vertices[grain][vertex][0]
//------------------------------------------------------------------------------
void grains::internal::write_grain_vertices(int id, double dx, double dy, std::ofstream& ofile, std::vector< std::vector <double> >& vertices){

   // only write on root process
   if(!vmpi::master) return;

   // write number of vertices
   ofile << vertices.size() << "\t";

   const int num_vertices = vertices.size();

   // write vertices
   for(int i = 0; i < num_vertices; i++){
      ofile << vertices[i][0]+dx << "\t" << vertices[i][1]+dy << "\t";
   }
   ofile << std::endl;

	return;
}
