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
#include "create.hpp"
#include "vmpi.hpp"
#include "vio.hpp"

// create module headers
#include "internal.hpp"

//------------------------------------------------------------------------------
// Function to write grain vertices to a file in the format
// num_vertices vx1 vy1 vx2 vy2 ...
//------------------------------------------------------------------------------
// Expects the vertices array in the format vertices[grain][vertex][0]
//------------------------------------------------------------------------------
void create::internal::write_grain_vertices(int id, double dx, double dy, std::ofstream& ofile, std::vector< std::vector <double> >& vertices){

   // only write on root process
   if(!vmpi::master) return;

   // write number of vertices
   ofile << vertices.size() << "\t";

   // write vertices
   for(int i = 0; i < vertices.size(); i++){
      ofile << vertices[i][0]+dx << "\t" << vertices[i][1]+dy << "\t";
   }
   ofile << std::endl;

   //----------------------------------------------------------------
   // commented out code to generate individual files for each grain
   //----------------------------------------------------------------
   //std::ofstream gfile;
   //std::stringstream fn;
   //fn << "grain" << id << ".txt";
   //gfile.open(fn.str());
   //for(int i = 0; i < vertices.size(); i++){
   //   gfile << vertices[i][0]+dx << "\t" << vertices[i][1]+dy << "\n";
   //}
   //gfile.close();

	return;
}
