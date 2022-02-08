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

//------------------------------------------------------------------------------
// Generates an idealised spinel structure
//
//------------------------------------------------------------------------------
void build_spinel(unitcell::unit_cell_t& unit_cell){

   // Set basic unit cell properties
   unit_cell.dimensions[0] = 1.0;
	unit_cell.dimensions[1] = 1.0;
	unit_cell.dimensions[2] = 1.0;

	unit_cell.shape[0][0] = 1;
	unit_cell.shape[0][1] = 0;
	unit_cell.shape[0][2] = 0;

	unit_cell.shape[1][0] = 0;
	unit_cell.shape[1][1] = 1;
	unit_cell.shape[1][2] = 0;

	unit_cell.shape[2][0] = 0;
	unit_cell.shape[2][1] = 0;
	unit_cell.shape[2][2] = 1;

	unit_cell.lcsize = 56;
	unit_cell.hcsize = 8;
	unit_cell.interaction_range = 1;
	unit_cell.atom.resize(56);
	unit_cell.surface_threshold = 8;

	//-----------------------------
	unit_cell.atom[0].x   = 0;
	unit_cell.atom[0].y   = 0;
	unit_cell.atom[0].z   = 0;
	unit_cell.atom[0].mat = 0;
	unit_cell.atom[0].lc  = 0;
	unit_cell.atom[0].hc  = 0;
	unit_cell.atom[0].ni  = 16;
   unit_cell.atom[0].nm  = false;
	//-----------------------------
	unit_cell.atom[1].x   = 0.75;
	unit_cell.atom[1].y   = 0.25;
	unit_cell.atom[1].z   = 0.25;
	unit_cell.atom[1].mat = 0;
	unit_cell.atom[1].lc  = 1;
	unit_cell.atom[1].hc  = 2;
	unit_cell.atom[1].ni  = 16;
   unit_cell.atom[1].nm  = false;
	//-----------------------------
	unit_cell.atom[2].x   = 0;
	unit_cell.atom[2].y   = 0.5;
	unit_cell.atom[2].z   = 0.5;
	unit_cell.atom[2].mat = 0;
	unit_cell.atom[2].lc  = 2;
	unit_cell.atom[2].hc  = 4;
	unit_cell.atom[2].ni  = 16;
   unit_cell.atom[2].nm  = false;
	//-----------------------------
	unit_cell.atom[3].x   = 0.75;
	unit_cell.atom[3].y   = 0.75;
	unit_cell.atom[3].z   = 0.75;
	unit_cell.atom[3].mat = 0;
	unit_cell.atom[3].lc  = 3;
	unit_cell.atom[3].hc  = 6;
	unit_cell.atom[3].ni  = 16;
   unit_cell.atom[3].nm  = false;
	//-----------------------------
	unit_cell.atom[4].x   = 0.5;
	unit_cell.atom[4].y   = 0;
	unit_cell.atom[4].z   = 0.5;
	unit_cell.atom[4].mat = 0;
	unit_cell.atom[4].lc  = 4;
	unit_cell.atom[4].hc  = 4;
	unit_cell.atom[4].ni  = 16;
   unit_cell.atom[4].nm  = false;
	//-----------------------------
	unit_cell.atom[5].x   = 0.25;
	unit_cell.atom[5].y   = 0.25;
	unit_cell.atom[5].z   = 0.75;
	unit_cell.atom[5].mat = 0;
	unit_cell.atom[5].lc  = 5;
	unit_cell.atom[5].hc  = 6;
	unit_cell.atom[5].ni  = 16;
   unit_cell.atom[5].nm  = false;
	//-----------------------------
	unit_cell.atom[6].x   = 0.5;
	unit_cell.atom[6].y   = 0.5;
	unit_cell.atom[6].z   = 0;
	unit_cell.atom[6].mat = 0;
	unit_cell.atom[6].lc  = 6;
	unit_cell.atom[6].hc  = 0;
	unit_cell.atom[6].ni  = 16;
   unit_cell.atom[6].nm  = false;
   //-----------------------------
	unit_cell.atom[7].x   = 0.25;
	unit_cell.atom[7].y   = 0.75;
	unit_cell.atom[7].z   = 0.25;
	unit_cell.atom[7].mat = 0;
	unit_cell.atom[7].lc  = 7;
	unit_cell.atom[7].hc  = 2;
	unit_cell.atom[7].ni  = 16;
   unit_cell.atom[7].nm  = false;
	//-----------------------------
	unit_cell.atom[8].x   = 0.375;
	unit_cell.atom[8].y   = 0.375;
	unit_cell.atom[8].z   = 0.375;
	unit_cell.atom[8].mat = 1;
	unit_cell.atom[8].lc  = 8;
	unit_cell.atom[8].hc  = 3;
	unit_cell.atom[8].ni  = 12;
   unit_cell.atom[8].nm  = false;
	//-----------------------------
	unit_cell.atom[9].x   = 0.125;
	unit_cell.atom[9].y   = 0.625;
	unit_cell.atom[9].z   = 0.875;
	unit_cell.atom[9].mat = 1;
	unit_cell.atom[9].lc  = 9;
	unit_cell.atom[9].hc  = 7;
	unit_cell.atom[9].ni  = 12;
   unit_cell.atom[9].nm  = false;
	//-----------------------------
	unit_cell.atom[10].x   = 0.625;
	unit_cell.atom[10].y   = 0.875;
	unit_cell.atom[10].z   = 0.125;
	unit_cell.atom[10].mat = 1;
	unit_cell.atom[10].lc  = 10;
	unit_cell.atom[10].hc  = 1;
	unit_cell.atom[10].ni  = 12;
   unit_cell.atom[10].nm  = false;
	//-----------------------------
	unit_cell.atom[11].x   = 0.875;
	unit_cell.atom[11].y   = 0.125;
	unit_cell.atom[11].z   = 0.625;
	unit_cell.atom[11].mat = 1;
	unit_cell.atom[11].lc  = 11;
	unit_cell.atom[11].hc  = 5;
	unit_cell.atom[11].ni  = 12;
   unit_cell.atom[11].nm  = false;
	//-----------------------------
	unit_cell.atom[12].x   = 0.375;
	unit_cell.atom[12].y   = 0.875;
	unit_cell.atom[12].z   = 0.875;
	unit_cell.atom[12].mat = 1;
	unit_cell.atom[12].lc  = 12;
	unit_cell.atom[12].hc  = 7;
	unit_cell.atom[12].ni  = 12;
   unit_cell.atom[12].nm  = false;
	//-----------------------------
	unit_cell.atom[13].x   = 0.125;
	unit_cell.atom[13].y   = 0.125;
	unit_cell.atom[13].z   = 0.375;
	unit_cell.atom[13].mat = 1;
	unit_cell.atom[13].lc  = 13;
	unit_cell.atom[13].hc  = 3;
	unit_cell.atom[13].ni  = 12;
   unit_cell.atom[13].nm  = false;
	//-----------------------------
	unit_cell.atom[14].x   = 0.625;
	unit_cell.atom[14].y   = 0.375;
	unit_cell.atom[14].z   = 0.625;
	unit_cell.atom[14].mat = 1;
	unit_cell.atom[14].lc  = 14;
	unit_cell.atom[14].hc  = 5;
	unit_cell.atom[14].ni  = 12;
   unit_cell.atom[14].nm  = false;
	//-----------------------------
	unit_cell.atom[15].x   = 0.875;
	unit_cell.atom[15].y   = 0.625;
	unit_cell.atom[15].z   = 0.125;
	unit_cell.atom[15].mat = 1;
	unit_cell.atom[15].lc  = 15;
	unit_cell.atom[15].hc  = 1;
	unit_cell.atom[15].ni  = 12;
   unit_cell.atom[15].nm  = false;
	//-----------------------------
	unit_cell.atom[16].x   = 0.875;
	unit_cell.atom[16].y   = 0.375;
	unit_cell.atom[16].z   = 0.875;
	unit_cell.atom[16].mat = 1;
	unit_cell.atom[16].lc  = 16;
	unit_cell.atom[16].hc  = 7;
	unit_cell.atom[16].ni  = 12;
   unit_cell.atom[16].nm  = false;
	//-----------------------------
	unit_cell.atom[17].x   = 0.625;
	unit_cell.atom[17].y   = 0.625;
	unit_cell.atom[17].z   = 0.375;
	unit_cell.atom[17].mat = 1;
	unit_cell.atom[17].lc  = 17;
	unit_cell.atom[17].hc  = 3;
	unit_cell.atom[17].ni  = 12;
   unit_cell.atom[17].nm  = false;
	//-----------------------------
	unit_cell.atom[18].x   = 0.125;
	unit_cell.atom[18].y   = 0.875;
	unit_cell.atom[18].z   = 0.625;
	unit_cell.atom[18].mat = 1;
	unit_cell.atom[18].lc  = 18;
	unit_cell.atom[18].hc  = 5;
	unit_cell.atom[18].ni  = 12;
   unit_cell.atom[18].nm  = false;
	//-----------------------------
	unit_cell.atom[19].x   = 0.375;
	unit_cell.atom[19].y   = 0.125;
	unit_cell.atom[19].z   = 0.125;
	unit_cell.atom[19].mat = 1;
	unit_cell.atom[19].lc  = 19;
	unit_cell.atom[19].hc  = 1;
	unit_cell.atom[19].ni  = 12;
   unit_cell.atom[19].nm  = false;
	//-----------------------------
	unit_cell.atom[20].x   = 0.875;
	unit_cell.atom[20].y   = 0.875;
	unit_cell.atom[20].z   = 0.375;
	unit_cell.atom[20].mat = 1;
	unit_cell.atom[20].lc  = 20;
	unit_cell.atom[20].hc  = 3;
	unit_cell.atom[20].ni  = 12;
   unit_cell.atom[20].nm  = false;
	//-----------------------------
	unit_cell.atom[21].x   = 0.625;
	unit_cell.atom[21].y   = 0.125;
	unit_cell.atom[21].z   = 0.875;
	unit_cell.atom[21].mat = 1;
	unit_cell.atom[21].lc  = 21;
	unit_cell.atom[21].hc  = 7;
	unit_cell.atom[21].ni  = 12;
   unit_cell.atom[21].nm  = false;
	//-----------------------------
	unit_cell.atom[22].x   = 0.125;
	unit_cell.atom[22].y   = 0.375;
	unit_cell.atom[22].z   = 0.125;
	unit_cell.atom[22].mat = 1;
	unit_cell.atom[22].lc  = 22;
	unit_cell.atom[22].hc  = 1;
	unit_cell.atom[22].ni  = 12;
   unit_cell.atom[22].nm  = false;
	//-----------------------------
	unit_cell.atom[23].x   = 0.375;
	unit_cell.atom[23].y   = 0.625;
	unit_cell.atom[23].z   = 0.625;
	unit_cell.atom[23].mat = 1;
	unit_cell.atom[23].lc  = 23;
	unit_cell.atom[23].hc  = 5;
	unit_cell.atom[23].ni  = 12;
   unit_cell.atom[23].nm  = false;
	//-----------------------------
	unit_cell.atom[24].x   = 0.125;
	unit_cell.atom[24].y   = 0.125;
	unit_cell.atom[24].z   = 0.125;
	unit_cell.atom[24].mat = 2;
	unit_cell.atom[24].lc  = 24;
	unit_cell.atom[24].hc  = 1;
	unit_cell.atom[24].ni  = 12;
   unit_cell.atom[24].nm  = false;
	//-----------------------------
	unit_cell.atom[25].x   = 0.375;
	unit_cell.atom[25].y   = 0.875;
	unit_cell.atom[25].z   = 0.625;
	unit_cell.atom[25].mat = 2;
	unit_cell.atom[25].lc  = 25;
	unit_cell.atom[25].hc  = 5;
	unit_cell.atom[25].ni  = 12;
   unit_cell.atom[25].nm  = false;
	//-----------------------------
	unit_cell.atom[26].x   = 0.875;
	unit_cell.atom[26].y   = 0.625;
	unit_cell.atom[26].z   = 0.375;
	unit_cell.atom[26].mat = 2;
	unit_cell.atom[26].lc  = 26;
	unit_cell.atom[26].hc  = 3;
	unit_cell.atom[26].ni  = 12;
   unit_cell.atom[26].nm  = false;
	//-----------------------------
	unit_cell.atom[27].x   = 0.625;
	unit_cell.atom[27].y   = 0.375;
	unit_cell.atom[27].z   = 0.875;
	unit_cell.atom[27].mat = 2;
	unit_cell.atom[27].lc  = 27;
	unit_cell.atom[27].hc  = 7;
	unit_cell.atom[27].ni  = 12;
   unit_cell.atom[27].nm  = false;
	//-----------------------------
	unit_cell.atom[28].x   = 0.875;
	unit_cell.atom[28].y   = 0.375;
	unit_cell.atom[28].z   = 0.125;
	unit_cell.atom[28].mat = 2;
	unit_cell.atom[28].lc  = 28;
	unit_cell.atom[28].hc  = 1;
	unit_cell.atom[28].ni  = 12;
   unit_cell.atom[28].nm  = false;
	//-----------------------------
	unit_cell.atom[29].x   = 0.625;
	unit_cell.atom[29].y   = 0.625;
	unit_cell.atom[29].z   = 0.625;
	unit_cell.atom[29].mat = 2;
	unit_cell.atom[29].lc  = 29;
	unit_cell.atom[29].hc  = 5;
	unit_cell.atom[29].ni  = 12;
   unit_cell.atom[29].nm  = false;
	//-----------------------------
	unit_cell.atom[30].x   = 0.375;
	unit_cell.atom[30].y   = 0.125;
	unit_cell.atom[30].z   = 0.875;
	unit_cell.atom[30].mat = 2;
	unit_cell.atom[30].lc  = 30;
	unit_cell.atom[30].hc  = 7;
	unit_cell.atom[30].ni  = 12;
   unit_cell.atom[30].nm  = false;
	//-----------------------------
	unit_cell.atom[31].x   = 0.125;
	unit_cell.atom[31].y   = 0.875;
	unit_cell.atom[31].z   = 0.375;
	unit_cell.atom[31].mat = 2;
	unit_cell.atom[31].lc  = 31;
	unit_cell.atom[31].hc  = 3;
	unit_cell.atom[31].ni  = 12;
   unit_cell.atom[31].nm  = false;
	//-----------------------------
	unit_cell.atom[32].x   = 0.125;
	unit_cell.atom[32].y   = 0.625;
	unit_cell.atom[32].z   = 0.625;
	unit_cell.atom[32].mat = 2;
	unit_cell.atom[32].lc  = 32;
	unit_cell.atom[32].hc  = 5;
	unit_cell.atom[32].ni  = 12;
   unit_cell.atom[32].nm  = false;
	//-----------------------------
	unit_cell.atom[33].x   = 0.375;
	unit_cell.atom[33].y   = 0.375;
	unit_cell.atom[33].z   = 0.125;
	unit_cell.atom[33].mat = 2;
	unit_cell.atom[33].lc  = 33;
	unit_cell.atom[33].hc  = 1;
	unit_cell.atom[33].ni  = 12;
   unit_cell.atom[33].nm  = false;
	//-----------------------------
	unit_cell.atom[34].x   = 0.875;
	unit_cell.atom[34].y   = 0.125;
	unit_cell.atom[34].z   = 0.875;
	unit_cell.atom[34].mat = 2;
	unit_cell.atom[34].lc  = 34;
	unit_cell.atom[34].hc  = 7;
	unit_cell.atom[34].ni  = 12;
   unit_cell.atom[34].nm  = false;
	//-----------------------------
	unit_cell.atom[35].x   = 0.625;
	unit_cell.atom[35].y   = 0.875;
	unit_cell.atom[35].z   = 0.375;
	unit_cell.atom[35].mat = 2;
	unit_cell.atom[35].lc  = 35;
	unit_cell.atom[35].hc  = 3;
	unit_cell.atom[35].ni  = 12;
   unit_cell.atom[35].nm  = false;
	//-----------------------------
	unit_cell.atom[36].x   = 0.875;
	unit_cell.atom[36].y   = 0.875;
	unit_cell.atom[36].z   = 0.625;
	unit_cell.atom[36].mat = 2;
	unit_cell.atom[36].lc  = 36;
	unit_cell.atom[36].hc  = 5;
	unit_cell.atom[36].ni  = 12;
   unit_cell.atom[36].nm  = false;
	//-----------------------------
	unit_cell.atom[37].x   = 0.625;
	unit_cell.atom[37].y   = 0.125;
	unit_cell.atom[37].z   = 0.125;
	unit_cell.atom[37].mat = 2;
	unit_cell.atom[37].lc  = 37;
	unit_cell.atom[37].hc  = 1;
	unit_cell.atom[37].ni  = 12;
   unit_cell.atom[37].nm  = false;
	//-----------------------------
	unit_cell.atom[38].x   = 0.375;
	unit_cell.atom[38].y   = 0.625;
	unit_cell.atom[38].z   = 0.375;
	unit_cell.atom[38].mat = 2;
	unit_cell.atom[38].lc  = 38;
	unit_cell.atom[38].hc  = 3;
	unit_cell.atom[38].ni  = 12;
   unit_cell.atom[38].nm  = false;
	//-----------------------------
	unit_cell.atom[39].x   = 0.125;
	unit_cell.atom[39].y   = 0.375;
	unit_cell.atom[39].z   = 0.875;
	unit_cell.atom[39].mat = 2;
	unit_cell.atom[39].lc  = 39;
	unit_cell.atom[39].hc  = 7;
	unit_cell.atom[39].ni  = 12;
   unit_cell.atom[39].nm  = false;
	//-----------------------------
	unit_cell.atom[40].x   = 0.625;
	unit_cell.atom[40].y   = 0.125;
	unit_cell.atom[40].z   = 0.625;
	unit_cell.atom[40].mat = 2;
	unit_cell.atom[40].lc  = 40;
	unit_cell.atom[40].hc  = 5;
	unit_cell.atom[40].ni  = 12;
   unit_cell.atom[40].nm  = false;
	//-----------------------------
	unit_cell.atom[41].x   = 0.875;
	unit_cell.atom[41].y   = 0.875;
	unit_cell.atom[41].z   = 0.125;
	unit_cell.atom[41].mat = 2;
	unit_cell.atom[41].lc  = 41;
	unit_cell.atom[41].hc  = 1;
	unit_cell.atom[41].ni  = 12;
   unit_cell.atom[41].nm  = false;
	//-----------------------------
	unit_cell.atom[42].x   = 0.375;
	unit_cell.atom[42].y   = 0.625;
	unit_cell.atom[42].z   = 0.875;
	unit_cell.atom[42].mat = 2;
	unit_cell.atom[42].lc  = 42;
	unit_cell.atom[42].hc  = 7;
	unit_cell.atom[42].ni  = 12;
   unit_cell.atom[42].nm  = false;
	//-----------------------------
	unit_cell.atom[43].x   = 0.125;
	unit_cell.atom[43].y   = 0.375;
	unit_cell.atom[43].z   = 0.375;
	unit_cell.atom[43].mat = 2;
	unit_cell.atom[43].lc  = 43;
	unit_cell.atom[43].hc  = 3;
	unit_cell.atom[43].ni  = 12;
   unit_cell.atom[43].nm  = false;
	//-----------------------------
	unit_cell.atom[44].x   = 0.375;
	unit_cell.atom[44].y   = 0.375;
	unit_cell.atom[44].z   = 0.625;
	unit_cell.atom[44].mat = 2;
	unit_cell.atom[44].lc  = 44;
	unit_cell.atom[44].hc  = 5;
	unit_cell.atom[44].ni  = 12;
   unit_cell.atom[44].nm  = false;
	//-----------------------------
	unit_cell.atom[45].x   = 0.125;
	unit_cell.atom[45].y   = 0.625;
	unit_cell.atom[45].z   = 0.125;
	unit_cell.atom[45].mat = 2;
	unit_cell.atom[45].lc  = 45;
	unit_cell.atom[45].hc  = 1;
	unit_cell.atom[45].ni  = 12;
   unit_cell.atom[45].nm  = false;
	//-----------------------------
	unit_cell.atom[46].x   = 0.875;
	unit_cell.atom[46].y   = 0.125;
	unit_cell.atom[46].z   = 0.375;
	unit_cell.atom[46].mat = 2;
	unit_cell.atom[46].lc  = 46;
	unit_cell.atom[46].hc  = 3;
	unit_cell.atom[46].ni  = 12;
   unit_cell.atom[46].nm  = false;
	//-----------------------------
	unit_cell.atom[47].x   = 0.625;
	unit_cell.atom[47].y   = 0.875;
	unit_cell.atom[47].z   = 0.875;
	unit_cell.atom[47].mat = 2;
	unit_cell.atom[47].lc  = 47;
	unit_cell.atom[47].hc  = 7;
	unit_cell.atom[47].ni  = 12;
   unit_cell.atom[47].nm  = false;
	//-----------------------------
	unit_cell.atom[48].x   = 0.625;
	unit_cell.atom[48].y   = 0.625;
	unit_cell.atom[48].z   = 0.125;
	unit_cell.atom[48].mat = 2;
	unit_cell.atom[48].lc  = 48;
	unit_cell.atom[48].hc  = 1;
	unit_cell.atom[48].ni  = 12;
   unit_cell.atom[48].nm  = false;
	//-----------------------------
	unit_cell.atom[49].x   = 0.875;
	unit_cell.atom[49].y   = 0.375;
	unit_cell.atom[49].z   = 0.625;
	unit_cell.atom[49].mat = 2;
	unit_cell.atom[49].lc  = 49;
	unit_cell.atom[49].hc  = 5;
	unit_cell.atom[49].ni  = 12;
   unit_cell.atom[49].nm  = false;
	//-----------------------------
	unit_cell.atom[50].x   = 0.375;
	unit_cell.atom[50].y   = 0.125;
	unit_cell.atom[50].z   = 0.375;
	unit_cell.atom[50].mat = 2;
	unit_cell.atom[50].lc  = 50;
	unit_cell.atom[50].hc  = 3;
	unit_cell.atom[50].ni  = 12;
   unit_cell.atom[50].nm  = false;
	//-----------------------------
	unit_cell.atom[51].x   = 0.125;
	unit_cell.atom[51].y   = 0.875;
	unit_cell.atom[51].z   = 0.875;
	unit_cell.atom[51].mat = 2;
	unit_cell.atom[51].lc  = 51;
	unit_cell.atom[51].hc  = 7;
	unit_cell.atom[51].ni  = 12;
   unit_cell.atom[51].nm  = false;
	//-----------------------------
	unit_cell.atom[52].x   = 0.375;
	unit_cell.atom[52].y   = 0.875;
	unit_cell.atom[52].z   = 0.125;
	unit_cell.atom[52].mat = 2;
	unit_cell.atom[52].lc  = 52;
	unit_cell.atom[52].hc  = 1;
	unit_cell.atom[52].ni  = 12;
   unit_cell.atom[52].nm  = false;
	//-----------------------------
	unit_cell.atom[53].x   = 0.125;
	unit_cell.atom[53].y   = 0.125;
	unit_cell.atom[53].z   = 0.625;
	unit_cell.atom[53].mat = 2;
	unit_cell.atom[53].lc  = 53;
	unit_cell.atom[53].hc  = 5;
	unit_cell.atom[53].ni  = 12;
   unit_cell.atom[53].nm  = false;
	//-----------------------------
	unit_cell.atom[54].x   = 0.875;
	unit_cell.atom[54].y   = 0.625;
	unit_cell.atom[54].z   = 0.875;
	unit_cell.atom[54].mat = 2;
	unit_cell.atom[54].lc  = 54;
	unit_cell.atom[54].hc  = 7;
	unit_cell.atom[54].ni  = 12;
   unit_cell.atom[54].nm  = false;
	//-----------------------------
	unit_cell.atom[55].x   = 0.625;
	unit_cell.atom[55].y   = 0.375;
	unit_cell.atom[55].z   = 0.375;
	unit_cell.atom[55].mat = 2;
	unit_cell.atom[55].lc  = 55;
	unit_cell.atom[55].hc  = 3;
	unit_cell.atom[55].ni  = 12;
   unit_cell.atom[55].nm  = false;

   unit_cell.cutoff_radius = 0.434; // normalised to unit cell size

   uc::internal::calculate_interactions(unit_cell);

   // Set actual unit cell size after calculating interactions
   unit_cell.dimensions[0] *= unitcell::internal::unit_cell_size_x;
   unit_cell.dimensions[1] *= unitcell::internal::unit_cell_size_y;
   unit_cell.dimensions[2] *= unitcell::internal::unit_cell_size_z;

   return;

}

} // end of internal namespace
} // end of unitcell namespace
