//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2019. All rights reserved.
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

void build_heusler(unitcell::unit_cell_t& unit_cell){

   // Set basic unit cell properties
   unit_cell.dimensions[0] = 1.0;
   unit_cell.dimensions[1] = 1.0;
   unit_cell.dimensions[2] = 1.0;

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
   unit_cell.hcsize=4;
   unit_cell.interaction_range=1;
   unit_cell.atom.resize(16);
   unit_cell.surface_threshold=8;
   //-----------------------------
   unit_cell.atom[0].x=0.0; // 4a
   unit_cell.atom[0].y=0.0;
   unit_cell.atom[0].z=0.0;
   unit_cell.atom[0].mat=0;
   unit_cell.atom[0].lc=0;
   unit_cell.atom[0].hc=0;
   unit_cell.atom[0].ni=26;
   //-----------------------------
   unit_cell.atom[1].x=0.5;
   unit_cell.atom[1].y=0.5;
   unit_cell.atom[1].z=0.0;
   unit_cell.atom[1].mat=0;
   unit_cell.atom[1].lc=0;
   unit_cell.atom[1].hc=0;
   unit_cell.atom[1].ni=26;
   //-----------------------------
   unit_cell.atom[2].x=0.5;
   unit_cell.atom[2].y=0.0;
   unit_cell.atom[2].z=0.5;
   unit_cell.atom[2].mat=0;
   unit_cell.atom[2].lc=0;
   unit_cell.atom[2].hc=1;
   unit_cell.atom[2].ni=26;
   //-----------------------------
   unit_cell.atom[3].x=0.0;
   unit_cell.atom[3].y=0.5;
   unit_cell.atom[3].z=0.5;
   unit_cell.atom[3].mat=0;
   unit_cell.atom[3].lc=0;
   unit_cell.atom[3].hc=1;
   unit_cell.atom[3].ni=26;
   //-----------------------------
   unit_cell.atom[4].x=0.5; // 4b
   unit_cell.atom[4].y=0.0;
   unit_cell.atom[4].z=0.0;
   unit_cell.atom[4].mat=1;
   unit_cell.atom[4].lc=1;
   unit_cell.atom[4].hc=0;
   unit_cell.atom[4].ni=26;
   //-----------------------------
   unit_cell.atom[5].x=0.0;
   unit_cell.atom[5].y=0.5;
   unit_cell.atom[5].z=0.0;
   unit_cell.atom[5].mat=1;
   unit_cell.atom[5].lc=1;
   unit_cell.atom[5].hc=0;
   unit_cell.atom[5].ni=26;
   //-----------------------------
   unit_cell.atom[6].x=0.0;
   unit_cell.atom[6].y=0.0;
   unit_cell.atom[6].z=0.5;
   unit_cell.atom[6].mat=1;
   unit_cell.atom[6].lc=1;
   unit_cell.atom[6].hc=1;
   unit_cell.atom[6].ni=26;
   //-----------------------------
   unit_cell.atom[7].x=0.5;
   unit_cell.atom[7].y=0.5;
   unit_cell.atom[7].z=0.5;
   unit_cell.atom[7].mat=1;
   unit_cell.atom[7].lc=1;
   unit_cell.atom[7].hc=1;
   unit_cell.atom[7].ni=26;
   //-----------------------------
   unit_cell.atom[8].x=0.25; // 4c
   unit_cell.atom[8].y=0.25;
   unit_cell.atom[8].z=0.25;
   unit_cell.atom[8].mat=2;
   unit_cell.atom[8].lc=2;
   unit_cell.atom[8].hc=2;
   unit_cell.atom[8].ni=26;
   //-----------------------------
   unit_cell.atom[9].x=0.75;
   unit_cell.atom[9].y=0.75;
   unit_cell.atom[9].z=0.25;
   unit_cell.atom[9].mat=2;
   unit_cell.atom[9].lc=2;
   unit_cell.atom[9].hc=2;
   unit_cell.atom[9].ni=26;
   //-----------------------------
   unit_cell.atom[10].x=0.75;
   unit_cell.atom[10].y=0.25;
   unit_cell.atom[10].z=0.75;
   unit_cell.atom[10].mat=2;
   unit_cell.atom[10].lc=2;
   unit_cell.atom[10].hc=3;
   unit_cell.atom[10].ni=26;
   //-----------------------------
   unit_cell.atom[11].x=0.25;
   unit_cell.atom[11].y=0.75;
   unit_cell.atom[11].z=0.75;
   unit_cell.atom[11].mat=2;
   unit_cell.atom[11].lc=2;
   unit_cell.atom[11].hc=3;
   unit_cell.atom[11].ni=26;
   //-----------------------------
   unit_cell.atom[12].x=0.75; // 4d
   unit_cell.atom[12].y=0.25;
   unit_cell.atom[12].z=0.25;
   unit_cell.atom[12].mat=3;
   unit_cell.atom[12].lc=3;
   unit_cell.atom[12].hc=2;
   unit_cell.atom[12].ni=26;
   //-----------------------------
   unit_cell.atom[13].x=0.25;
   unit_cell.atom[13].y=0.75;
   unit_cell.atom[13].z=0.25;
   unit_cell.atom[13].mat=3;
   unit_cell.atom[13].lc=3;
   unit_cell.atom[13].hc=2;
   unit_cell.atom[13].ni=26;
   //-----------------------------
   unit_cell.atom[14].x=0.25;
   unit_cell.atom[14].y=0.25;
   unit_cell.atom[14].z=0.75;
   unit_cell.atom[14].mat=3;
   unit_cell.atom[14].lc=3;
   unit_cell.atom[14].hc=3;
   unit_cell.atom[14].ni=26;
   //-----------------------------
   unit_cell.atom[15].x=0.75;
   unit_cell.atom[15].y=0.75;
   unit_cell.atom[15].z=0.75;
   unit_cell.atom[15].mat=3;
   unit_cell.atom[15].lc=3;
   unit_cell.atom[15].hc=3;
   unit_cell.atom[15].ni=26;

   unit_cell.cutoff_radius = 0.708; // normalised to unit cell size

   uc::internal::calculate_interactions(unit_cell);

   // Set actual unit cell size after calculating interactions
   unit_cell.dimensions[0] *= unitcell::internal::unit_cell_size_x;
   unit_cell.dimensions[1] *= unitcell::internal::unit_cell_size_y;
   unit_cell.dimensions[2] *= unitcell::internal::unit_cell_size_z;

   return;

}

} // end of internal namespace
} // end of unitcell namespace
