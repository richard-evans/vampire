//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Jack B Collings 2021. All rights reserved.
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

// number of interactions notional value used.

namespace unitcell{
namespace internal{



void build_SmFeN(unitcell::unit_cell_t& unit_cell){

   // Set basic unit cell properties
   unit_cell.dimensions[0] = 1.0;
   unit_cell.dimensions[1] = 1.0;
   unit_cell.dimensions[2] = 4.802 / 8.566;

   unit_cell.shape[0][0] = 1;
   unit_cell.shape[0][1] = 0;
   unit_cell.shape[0][2] = 0;

   unit_cell.shape[1][0] = 0;
   unit_cell.shape[1][1] = 1;
   unit_cell.shape[1][2] = 0;

   unit_cell.shape[2][0] = 0;
   unit_cell.shape[2][1] = 0;
   unit_cell.shape[2][2] = 1;

   unit_cell.lcsize = 28;
   unit_cell.hcsize = 6;            // 
   unit_cell.interaction_range = 1;
   unit_cell.atom.resize(28);
   unit_cell.surface_threshold = 10; // Not yet properly calculated. Also atom.ni not properly calculated, but ni is now redundant.

   //-----------------------------
   unit_cell.atom[0].x     = 0.0000;
   unit_cell.atom[0].y     = 0.0000;
   unit_cell.atom[0].z     = 0.0000;
   unit_cell.atom[0].mat   = 0;
   unit_cell.atom[0].lc    = 0;
   unit_cell.atom[0].hc    = 0;
   unit_cell.atom[0].nm    = false;
   unit_cell.atom[0].ni    = 10;
   //-----------------------------
   unit_cell.atom[1].x     = 0.5000;
   unit_cell.atom[1].y     = 0.5000;
   unit_cell.atom[1].z     = 0.5000;
   unit_cell.atom[1].mat   = 0;
   unit_cell.atom[1].lc    = 1;
   unit_cell.atom[1].hc    = 2;
   unit_cell.atom[1].nm    = false;
   unit_cell.atom[1].ni    = 10;
   //-----------------------------
   unit_cell.atom[2].x     = 0.3534;
   unit_cell.atom[2].y     = 0.0000;
   unit_cell.atom[2].z     = 0.0000;
   unit_cell.atom[2].mat   = 1;
   unit_cell.atom[2].lc    = 2;
   unit_cell.atom[2].hc    = 0;
   unit_cell.atom[2].nm    = false;
   unit_cell.atom[2].ni    = 10;
   //-----------------------------
   unit_cell.atom[3].x     = 0.6466;
   unit_cell.atom[3].y     = 0.0000;
   unit_cell.atom[3].z     = 0.0000;
   unit_cell.atom[3].mat   = 1;
   unit_cell.atom[3].lc    = 3;
   unit_cell.atom[3].hc    = 0;
   unit_cell.atom[3].nm    = false;
   unit_cell.atom[3].ni    = 10;
   //-----------------------------
   unit_cell.atom[4].x     = 0.0000;
   unit_cell.atom[4].y     = 0.3534;
   unit_cell.atom[4].z     = 0.0000;
   unit_cell.atom[4].mat   = 1;
   unit_cell.atom[4].lc    = 4;
   unit_cell.atom[4].hc    = 0;
   unit_cell.atom[4].nm    = false;
   unit_cell.atom[4].ni    = 10;
   //-----------------------------
   unit_cell.atom[5].x     = 0.0000;
   unit_cell.atom[5].y     = 0.6466;
   unit_cell.atom[5].z     = 0.0000;
   unit_cell.atom[5].mat   = 1;
   unit_cell.atom[5].lc    = 5;
   unit_cell.atom[5].hc    = 0;
   unit_cell.atom[5].nm    = false;
   unit_cell.atom[5].ni    = 10;
   //----------------------------
   unit_cell.atom[6].x     = 0.2753;
   unit_cell.atom[6].y     = 0.5000;
   unit_cell.atom[6].z     = 0.0000;
   unit_cell.atom[6].mat   = 1;
   unit_cell.atom[6].lc    = 6;
   unit_cell.atom[6].hc    = 0;
   unit_cell.atom[6].nm    = false;
   unit_cell.atom[6].ni    = 10;
   //----------------------------
   unit_cell.atom[7].x     = 0.7247;
   unit_cell.atom[7].y     = 0.5000;
   unit_cell.atom[7].z     = 0.0000;
   unit_cell.atom[7].mat   = 1;
   unit_cell.atom[7].lc    = 7;
   unit_cell.atom[7].hc    = 0;
   unit_cell.atom[7].nm    = false;
   unit_cell.atom[7].ni    = 10;
   //----------------------------
   unit_cell.atom[8].x     = 0.5000;
   unit_cell.atom[8].y     = 0.2753;
   unit_cell.atom[8].z     = 0.0000;
   unit_cell.atom[8].mat   = 1;
   unit_cell.atom[8].lc    = 8;
   unit_cell.atom[8].hc    = 0;
   unit_cell.atom[8].nm    = false;
   unit_cell.atom[8].ni    = 10;
   //----------------------------
   unit_cell.atom[9].x     = 0.5000;
   unit_cell.atom[9].y     = 0.7247;
   unit_cell.atom[9].z     = 0.0000;
   unit_cell.atom[9].mat   = 1;
   unit_cell.atom[9].lc    = 9;
   unit_cell.atom[9].hc    = 0;
   unit_cell.atom[9].nm    = false;
   unit_cell.atom[9].ni    = 10;
   //----------------------------
   unit_cell.atom[10].x    = 0.2500;
   unit_cell.atom[10].y    = 0.2500;
   unit_cell.atom[10].z    = 0.2500;
   unit_cell.atom[10].mat  = 1;
   unit_cell.atom[10].lc   = 10;
   unit_cell.atom[10].hc   = 1;
   unit_cell.atom[10].nm   = false;
   unit_cell.atom[10].ni   = 10;
   //----------------------------
   unit_cell.atom[11].x    = 0.7500;
   unit_cell.atom[11].y    = 0.7500;
   unit_cell.atom[11].z    = 0.2500;
   unit_cell.atom[11].mat  = 1;
   unit_cell.atom[11].lc   = 11;
   unit_cell.atom[11].hc   = 1;
   unit_cell.atom[11].nm   = false;
   unit_cell.atom[11].ni   = 10;
   //----------------------------
   unit_cell.atom[12].x    = 0.7500;
   unit_cell.atom[12].y    = 0.2500;
   unit_cell.atom[12].z    = 0.2500;
   unit_cell.atom[12].mat  = 1;
   unit_cell.atom[12].lc   = 12;
   unit_cell.atom[12].hc   = 1;
   unit_cell.atom[12].nm   = false;
   unit_cell.atom[12].ni   = 10;
   //----------------------------
   unit_cell.atom[13].x    = 0.2500;
   unit_cell.atom[13].y    = 0.7500;
   unit_cell.atom[13].z    = 0.2500;
   unit_cell.atom[13].mat  = 1;
   unit_cell.atom[13].lc   = 13;
   unit_cell.atom[13].hc   = 1;
   unit_cell.atom[13].nm   = false;
   unit_cell.atom[13].ni   = 10;
   //----------------------------
   unit_cell.atom[14].x    = 0.8534;
   unit_cell.atom[14].y    = 0.5000;
   unit_cell.atom[14].z    = 0.5000;
   unit_cell.atom[14].mat  = 1;
   unit_cell.atom[14].lc   = 14;
   unit_cell.atom[14].hc   = 2;
   unit_cell.atom[14].nm   = false;
   unit_cell.atom[14].ni   = 10;
   //----------------------------
   unit_cell.atom[15].x    = 0.1466;
   unit_cell.atom[15].y    = 0.5000;
   unit_cell.atom[15].z    = 0.5000;
   unit_cell.atom[15].mat  = 1;
   unit_cell.atom[15].lc   = 15;
   unit_cell.atom[15].hc   = 2;
   unit_cell.atom[15].nm   = false;
   unit_cell.atom[15].ni   = 10;
   //----------------------------
   unit_cell.atom[16].x    = 0.5000;
   unit_cell.atom[16].y    = 0.8534;
   unit_cell.atom[16].z    = 0.5000;
   unit_cell.atom[16].mat  = 1;
   unit_cell.atom[16].lc   = 16;
   unit_cell.atom[16].hc   = 2;
   unit_cell.atom[16].nm   = false;
   unit_cell.atom[16].ni   = 10;
   //----------------------------
   unit_cell.atom[17].x    = 0.5000;
   unit_cell.atom[17].y    = 0.1466;
   unit_cell.atom[17].z    = 0.5000;
   unit_cell.atom[17].mat  = 1;
   unit_cell.atom[17].lc   = 17;
   unit_cell.atom[17].hc   = 2;
   unit_cell.atom[17].nm   = false;
   unit_cell.atom[17].ni   = 10;
   //----------------------------
   unit_cell.atom[18].x    = 0.7753;
   unit_cell.atom[18].y    = 0.0000;
   unit_cell.atom[18].z    = 0.5000;
   unit_cell.atom[18].mat  = 1;
   unit_cell.atom[18].lc   = 18;
   unit_cell.atom[18].hc   = 2;
   unit_cell.atom[18].nm   = false;
   unit_cell.atom[18].ni   = 10;
   //----------------------------
   unit_cell.atom[19].x    = 0.2247;
   unit_cell.atom[19].y    = 0.0000;
   unit_cell.atom[19].z    = 0.5000;
   unit_cell.atom[19].mat  = 1;
   unit_cell.atom[19].lc   = 19;
   unit_cell.atom[19].hc   = 2;
   unit_cell.atom[19].nm   = false;
   unit_cell.atom[19].ni   = 10;
   //----------------------------
   unit_cell.atom[20].x    = 0.0000;
   unit_cell.atom[20].y    = 0.7753;
   unit_cell.atom[20].z    = 0.5000;
   unit_cell.atom[20].mat  = 1;
   unit_cell.atom[20].lc   = 20;
   unit_cell.atom[20].hc   = 2;
   unit_cell.atom[20].nm   = false;
   unit_cell.atom[20].ni   = 10;
   //----------------------------
   unit_cell.atom[21].x    = 0.0000;
   unit_cell.atom[21].y    = 0.2247;
   unit_cell.atom[21].z    = 0.5000;
   unit_cell.atom[21].mat  = 1;
   unit_cell.atom[21].lc   = 21;
   unit_cell.atom[21].hc   = 2;
   unit_cell.atom[21].nm   = false;
   unit_cell.atom[21].ni   = 10;
   //----------------------------
   unit_cell.atom[22].x    = 0.7500;
   unit_cell.atom[22].y    = 0.7500;
   unit_cell.atom[22].z    = 0.7500;
   unit_cell.atom[22].mat  = 1;
   unit_cell.atom[22].lc   = 22;
   unit_cell.atom[22].hc   = 3;
   unit_cell.atom[22].nm   = false;
   unit_cell.atom[22].ni   = 10;
   //----------------------------
   unit_cell.atom[23].x    = 0.2500;
   unit_cell.atom[23].y    = 0.2500;
   unit_cell.atom[23].z    = 0.7500;
   unit_cell.atom[23].mat  = 1;
   unit_cell.atom[23].lc   = 23;
   unit_cell.atom[23].hc   = 3;
   unit_cell.atom[23].nm   = false;
   unit_cell.atom[23].ni   = 10;
   //----------------------------
   unit_cell.atom[24].x    = 0.2500;
   unit_cell.atom[24].y    = 0.7500;
   unit_cell.atom[24].z    = 0.7500;
   unit_cell.atom[24].mat  = 1;
   unit_cell.atom[24].lc   = 24;
   unit_cell.atom[24].hc   = 3;
   unit_cell.atom[24].nm   = false;
   unit_cell.atom[24].ni   = 10;
   //----------------------------
   unit_cell.atom[25].x    = 0.7500;
   unit_cell.atom[25].y    = 0.2500;
   unit_cell.atom[25].z    = 0.7500;
   unit_cell.atom[25].mat  = 1;
   unit_cell.atom[25].lc   = 25;
   unit_cell.atom[25].hc   = 3;
   unit_cell.atom[25].nm   = false;
   unit_cell.atom[25].ni   = 10;
   //----------------------------
   unit_cell.atom[26].x    = 0.0000;
   unit_cell.atom[26].y    = 0.0000;
   unit_cell.atom[26].z    = 0.5000;
   unit_cell.atom[26].mat  = 2;
   unit_cell.atom[26].lc   = 26;
   unit_cell.atom[26].hc   = 2;
   unit_cell.atom[26].nm   = true; //Nitrogen
   unit_cell.atom[26].ni   = 10;
   //----------------------------
   unit_cell.atom[27].x    = 0.5000;
   unit_cell.atom[27].y    = 0.5000;
   unit_cell.atom[27].z    = 0.0000;
   unit_cell.atom[27].mat  = 2;
   unit_cell.atom[27].lc   = 27;
   unit_cell.atom[27].hc   = 0;
   unit_cell.atom[27].nm   = true;
   unit_cell.atom[27].ni   = 10;

   unit_cell.cutoff_radius = 0.50; // normalised to x-axis unit cell length

   uc::internal::calculate_interactions(unit_cell);

   // Set actual unit cell size after calculating interactions
   unit_cell.dimensions[0] *= unitcell::internal::unit_cell_size_x;
   unit_cell.dimensions[1] *= unitcell::internal::unit_cell_size_y;
   unit_cell.dimensions[2] *= unitcell::internal::unit_cell_size_z;

   return;

}

} // End internal namespace
} // End unitcell namespace
