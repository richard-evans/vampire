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
#include <cmath>

// Vampire headers
#include "unitcell.hpp"

// unitcell module headers
#include "internal.hpp"

namespace unitcell{
namespace internal{

void build_kagome(unitcell::unit_cell_t& unit_cell){

   // rescale y and z unit cell size to preserve unit cell dimensions
   unit_cell.dimensions[0]=1.0;
   unit_cell.dimensions[1]=1.0/3.0*sqrt(3.0);
   unit_cell.dimensions[2]=2.0*sqrt(2.0)/3.0;

   unit_cell.shape[0][0]=1.0;
   unit_cell.shape[0][1]=0.0;
   unit_cell.shape[0][2]=0.0;

   unit_cell.shape[1][0]=0.0;
   unit_cell.shape[1][1]=1.0;
   unit_cell.shape[1][2]=0.0;

   unit_cell.shape[2][0]=0.0;
   unit_cell.shape[2][1]=0.0;
   unit_cell.shape[2][2]=1.0;

   unit_cell.lcsize=4;
   unit_cell.hcsize=1;
   unit_cell.interaction_range=1;
   unit_cell.atom.resize(4);
   unit_cell.surface_threshold=3;
   //-----------------------------
   unit_cell.atom[0].x = 0.0;
   unit_cell.atom[0].y = 0.0;
   unit_cell.atom[0].z = 0.0;
   unit_cell.atom[0].lc = 0;
   unit_cell.atom[0].hc = 0;
   unit_cell.atom[0].ni = 3;
   //-----------------------------
   unit_cell.atom[1].x = 1.0/3.0;
   unit_cell.atom[1].y = 0.0;
   unit_cell.atom[1].z = 0.0;
   unit_cell.atom[1].lc = 1;
   unit_cell.atom[1].hc = 0;
   unit_cell.atom[1].ni = 3;
   //-----------------------------
   unit_cell.atom[2].x = 0.5;
   unit_cell.atom[2].y = 0.5;
   unit_cell.atom[2].z = 0.0;
   unit_cell.atom[2].lc = 0;
   unit_cell.atom[2].hc = 0;
   unit_cell.atom[2].ni = 3;
   //-----------------------------
   unit_cell.atom[3].x = 5.0/6.0;
   unit_cell.atom[3].y = 0.5;
   unit_cell.atom[3].z = 0.0;
   unit_cell.atom[3].lc = 3;
   unit_cell.atom[3].hc = 0;
   unit_cell.atom[3].ni = 3;

   unit_cell.cutoff_radius = 0.334; // normalised to unit cell size

   uc::internal::calculate_interactions(unit_cell);

   // Set actual unit cell size after calculating interactions
   unit_cell.dimensions[0] *= unitcell::internal::unit_cell_size_x;
   unit_cell.dimensions[1] *= unitcell::internal::unit_cell_size_y;
   unit_cell.dimensions[2] *= unitcell::internal::unit_cell_size_z;

   return;

}

} // end of internal namespace
} // end of unitcell namespace
