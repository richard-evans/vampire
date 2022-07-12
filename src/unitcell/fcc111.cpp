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

//------------------------------------------------------------------------------
// Generate a face-centred crystal with 111 orientation perpendicular to z
//------------------------------------------------------------------------------
void build_face_centred_cubic_111(unitcell::unit_cell_t& unit_cell){


   // Set basic unit cell properties
   unit_cell.dimensions[0] = 1.0;
	unit_cell.dimensions[1] = 1.0 / sqrt(3.0);
	unit_cell.dimensions[2] = sqrt(2.0)/2.0;

	unit_cell.shape[0][0] = 1;
	unit_cell.shape[0][1] = 0;
	unit_cell.shape[0][2] = 0;

	unit_cell.shape[1][0] = 0;
	unit_cell.shape[1][1] = 1;
	unit_cell.shape[1][2] = 0;

	unit_cell.shape[2][0] = 0;
	unit_cell.shape[2][1] = 0;
	unit_cell.shape[2][2] = 1;

	unit_cell.lcsize = 4;
	unit_cell.hcsize = 3;
	unit_cell.interaction_range = 1;
	unit_cell.atom.resize(24);
	unit_cell.surface_threshold = 12;

	//-----------------------------
	unit_cell.atom[0].x   = 0;
	unit_cell.atom[0].y   = 0;
	unit_cell.atom[0].z   = 0;
	unit_cell.atom[0].mat = uc::internal::sublattice_materials ? 3 : 0;
	unit_cell.atom[0].lc  = 3;
	unit_cell.atom[0].hc  = 0;
	unit_cell.atom[0].ni  = 12;
	unit_cell.atom[0].nm  = 0;
	//-----------------------------
	unit_cell.atom[1].x   = 0;
	unit_cell.atom[1].y   = 0.5;
	unit_cell.atom[1].z   = 0;
	unit_cell.atom[1].mat = uc::internal::sublattice_materials ? 2 : 0;
	unit_cell.atom[1].lc  = 2;
	unit_cell.atom[1].hc  = 0;
	unit_cell.atom[1].ni  = 12;
	unit_cell.atom[1].nm  = 0;
	//-----------------------------
	unit_cell.atom[2].x   = 0.0833333333;
	unit_cell.atom[2].y   = 0.25;
	unit_cell.atom[2].z   = 0.333333333;
	unit_cell.atom[2].mat = uc::internal::sublattice_materials ? 0 : 0;
	unit_cell.atom[2].lc  = 0;
	unit_cell.atom[2].hc  = 1;
	unit_cell.atom[2].ni  = 12;
	unit_cell.atom[2].nm  = 0;
	//-----------------------------
	unit_cell.atom[3].x   = 0.0833333333;
	unit_cell.atom[3].y   = 0.75;
	unit_cell.atom[3].z   = 0.333333333;
	unit_cell.atom[3].mat = uc::internal::sublattice_materials ? 1 : 0;
	unit_cell.atom[3].lc  = 1;
	unit_cell.atom[3].hc  = 1;
	unit_cell.atom[3].ni  = 12;
	unit_cell.atom[3].nm  = 0;
	//-----------------------------
	unit_cell.atom[4].x   = 0.166666667;
	unit_cell.atom[4].y   = 0;
	unit_cell.atom[4].z   = 0.666666667;
	unit_cell.atom[4].mat = uc::internal::sublattice_materials ? 2 : 0;
	unit_cell.atom[4].lc  = 2;
	unit_cell.atom[4].hc  = 2;
	unit_cell.atom[4].ni  = 12;
	unit_cell.atom[4].nm  = 0;
	//-----------------------------
	unit_cell.atom[5].x   = 0.166666667;
	unit_cell.atom[5].y   = 0.5;
	unit_cell.atom[5].z   = 0.666666667;
	unit_cell.atom[5].mat = uc::internal::sublattice_materials ? 3 : 0;
	unit_cell.atom[5].lc  = 3;
	unit_cell.atom[5].hc  = 2;
	unit_cell.atom[5].ni  = 12;
	unit_cell.atom[5].nm  = 0;
	//-----------------------------
	unit_cell.atom[6].x   = 0.25;
	unit_cell.atom[6].y   = 0.25;
	unit_cell.atom[6].z   = 0;
	unit_cell.atom[6].mat = uc::internal::sublattice_materials ? 1 : 0;
	unit_cell.atom[6].lc  = 1;
	unit_cell.atom[6].hc  = 0;
	unit_cell.atom[6].ni  = 12;
	unit_cell.atom[6].nm  = 0;
	//-----------------------------
	unit_cell.atom[7].x   = 0.25;
	unit_cell.atom[7].y   = 0.75;
	unit_cell.atom[7].z   = 0;
	unit_cell.atom[7].mat = uc::internal::sublattice_materials ? 0 : 0;
	unit_cell.atom[7].lc  = 0;
	unit_cell.atom[7].hc  = 0;
	unit_cell.atom[7].ni  = 12;
	unit_cell.atom[7].nm  = 0;
	//-----------------------------
	unit_cell.atom[8].x   = 0.333333333;
	unit_cell.atom[8].y   = 0;
	unit_cell.atom[8].z   = 0.333333333;
	unit_cell.atom[8].mat = uc::internal::sublattice_materials ? 3 : 0;
	unit_cell.atom[8].lc  = 3;
	unit_cell.atom[8].hc  = 1;
	unit_cell.atom[8].ni  = 12;
	unit_cell.atom[8].nm  = 0;
	//-----------------------------
	unit_cell.atom[9].x   = 0.333333333;
	unit_cell.atom[9].y   = 0.5;
	unit_cell.atom[9].z   = 0.333333333;
	unit_cell.atom[9].mat = uc::internal::sublattice_materials ? 2 : 0;
	unit_cell.atom[9].lc  = 2;
	unit_cell.atom[9].hc  = 1;
	unit_cell.atom[9].ni  = 12;
	unit_cell.atom[9].nm  = 0;
	//-----------------------------
	unit_cell.atom[10].x   = 0.416666667;
	unit_cell.atom[10].y   = 0.25;
	unit_cell.atom[10].z   = 0.666666667;
	unit_cell.atom[10].mat = uc::internal::sublattice_materials ? 0 : 0;
	unit_cell.atom[10].lc  = 0;
	unit_cell.atom[10].hc  = 2;
	unit_cell.atom[10].ni  = 12;
	unit_cell.atom[10].nm  = 0;
	//-----------------------------
	unit_cell.atom[11].x   = 0.416666667;
	unit_cell.atom[11].y   = 0.75;
	unit_cell.atom[11].z   = 0.666666667;
	unit_cell.atom[11].mat = uc::internal::sublattice_materials ? 1 : 0;
	unit_cell.atom[11].lc  = 1;
	unit_cell.atom[11].hc  = 2;
	unit_cell.atom[11].ni  = 12;
	unit_cell.atom[11].nm  = 0;
	//-----------------------------
	unit_cell.atom[12].x   = 0.5;
	unit_cell.atom[12].y   = 0;
	unit_cell.atom[12].z   = 0;
	unit_cell.atom[12].mat = uc::internal::sublattice_materials ? 2 : 0;
	unit_cell.atom[12].lc  = 2;
	unit_cell.atom[12].hc  = 0;
	unit_cell.atom[12].ni  = 12;
	unit_cell.atom[12].nm  = 0;
	//-----------------------------
	unit_cell.atom[13].x   = 0.5;
	unit_cell.atom[13].y   = 0.5;
	unit_cell.atom[13].z   = 0;
	unit_cell.atom[13].mat = uc::internal::sublattice_materials ? 3 : 0;
	unit_cell.atom[13].lc  = 3;
	unit_cell.atom[13].hc  = 0;
	unit_cell.atom[13].ni  = 12;
	unit_cell.atom[13].nm  = 0;
	//-----------------------------
	unit_cell.atom[14].x   = 0.583333333;
	unit_cell.atom[14].y   = 0.25;
	unit_cell.atom[14].z   = 0.333333333;
	unit_cell.atom[14].mat = uc::internal::sublattice_materials ? 1 : 0;
	unit_cell.atom[14].lc  = 1;
	unit_cell.atom[14].hc  = 1;
	unit_cell.atom[14].ni  = 12;
	unit_cell.atom[14].nm  = 0;
	//-----------------------------
	unit_cell.atom[15].x   = 0.583333333;
	unit_cell.atom[15].y   = 0.75;
	unit_cell.atom[15].z   = 0.333333333;
	unit_cell.atom[15].mat = uc::internal::sublattice_materials ? 0 : 0;
	unit_cell.atom[15].lc  = 0;
	unit_cell.atom[15].hc  = 1;
	unit_cell.atom[15].ni  = 12;
	unit_cell.atom[15].nm  = 0;
	//-----------------------------
	unit_cell.atom[16].x   = 0.666666667;
	unit_cell.atom[16].y   = 0;
	unit_cell.atom[16].z   = 0.666666667;
	unit_cell.atom[16].mat = uc::internal::sublattice_materials ? 3 : 0;
	unit_cell.atom[16].lc  = 3;
	unit_cell.atom[16].hc  = 2;
	unit_cell.atom[16].ni  = 12;
	unit_cell.atom[16].nm  = 0;
	//-----------------------------
	unit_cell.atom[17].x   = 0.666666667;
	unit_cell.atom[17].y   = 0.5;
	unit_cell.atom[17].z   = 0.666666667;
	unit_cell.atom[17].mat = uc::internal::sublattice_materials ? 2 : 0;
	unit_cell.atom[17].lc  = 2;
	unit_cell.atom[17].hc  = 2;
	unit_cell.atom[17].ni  = 12;
	unit_cell.atom[17].nm  = 0;
	//-----------------------------
	unit_cell.atom[18].x   = 0.75;
	unit_cell.atom[18].y   = 0.25;
	unit_cell.atom[18].z   = 0;
	unit_cell.atom[18].mat = uc::internal::sublattice_materials ? 0 : 0;
	unit_cell.atom[18].lc  = 0;
	unit_cell.atom[18].hc  = 0;
	unit_cell.atom[18].ni  = 12;
	unit_cell.atom[18].nm  = 0;
	//-----------------------------
	unit_cell.atom[19].x   = 0.75;
	unit_cell.atom[19].y   = 0.75;
	unit_cell.atom[19].z   = 0;
	unit_cell.atom[19].mat = uc::internal::sublattice_materials ? 1 : 0;
	unit_cell.atom[19].lc  = 1;
	unit_cell.atom[19].hc  = 0;
	unit_cell.atom[19].ni  = 12;
	unit_cell.atom[19].nm  = 0;
	//-----------------------------
	unit_cell.atom[20].x   = 0.833333333;
	unit_cell.atom[20].y   = 0;
	unit_cell.atom[20].z   = 0.333333333;
	unit_cell.atom[20].mat = uc::internal::sublattice_materials ? 2 : 0;
	unit_cell.atom[20].lc  = 2;
	unit_cell.atom[20].hc  = 1;
	unit_cell.atom[20].ni  = 12;
	unit_cell.atom[20].nm  = 0;
	//-----------------------------
	unit_cell.atom[21].x   = 0.833333333;
	unit_cell.atom[21].y   = 0.5;
	unit_cell.atom[21].z   = 0.333333333;
	unit_cell.atom[21].mat = uc::internal::sublattice_materials ? 3 : 0;
	unit_cell.atom[21].lc  = 3;
	unit_cell.atom[21].hc  = 1;
	unit_cell.atom[21].ni  = 12;
	unit_cell.atom[21].nm  = 0;
	//-----------------------------
	unit_cell.atom[22].x   = 0.916666667;
	unit_cell.atom[22].y   = 0.25;
	unit_cell.atom[22].z   = 0.666666667;
	unit_cell.atom[22].mat = uc::internal::sublattice_materials ? 1 : 0;
	unit_cell.atom[22].lc  = 1;
	unit_cell.atom[22].hc  = 2;
	unit_cell.atom[22].ni  = 12;
	unit_cell.atom[22].nm  = 0;
	//-----------------------------
	unit_cell.atom[23].x   = 0.916666667;
	unit_cell.atom[23].y   = 0.75;
	unit_cell.atom[23].z   = 0.666666667;
	unit_cell.atom[23].mat = uc::internal::sublattice_materials ? 0 : 0;
	unit_cell.atom[23].lc  = 0;
	unit_cell.atom[23].hc  = 2;
	unit_cell.atom[23].ni  = 12;
	unit_cell.atom[23].nm  = 0;

   unit_cell.cutoff_radius = sqrt(0.3333333333333)*0.5; // normalised to unit cell size

   uc::internal::calculate_interactions(unit_cell);

   // Set actual unit cell size after calculating interactions
   unit_cell.dimensions[0] *= unitcell::internal::unit_cell_size_x;
   unit_cell.dimensions[1] *= unitcell::internal::unit_cell_size_y;
   unit_cell.dimensions[2] *= unitcell::internal::unit_cell_size_z;

   return;

}

} // end of internal namespace
} // end of unitcell namespace
