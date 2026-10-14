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
#include <iostream>

// Vampire headers
#include "unitcell.hpp"

// unitcell module headers
#include "internal.hpp"

namespace unitcell{
namespace internal{

void build_rock_salt_supercell(unitcell::unit_cell_t& unit_cell){

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

   unit_cell.lcsize=3;
   unit_cell.hcsize=4;
   unit_cell.interaction_range=1;
   unit_cell.atom.resize(64);
   unit_cell.surface_threshold=6;
   //-----------------------------
   unit_cell.atom[0].x = 0;
   unit_cell.atom[0].y = 0;
   unit_cell.atom[0].z = 0;
   unit_cell.atom[0].mat = uc::internal::sublattice_materials ? 0 : 1; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[0].lc = 1;
   unit_cell.atom[0].hc = 0;
   unit_cell.atom[0].ni = 6;
   //---------------------------------------------
   unit_cell.atom[1].x = 0.25;
   unit_cell.atom[1].y = 0;
   unit_cell.atom[1].z = 0;
   unit_cell.atom[1].mat = uc::internal::sublattice_materials ? 1 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[1].lc = 2;
   unit_cell.atom[1].hc = 0;
   unit_cell.atom[1].ni = 6;
   //---------------------------------------------
   unit_cell.atom[2].x = 0;
   unit_cell.atom[2].y = 0.25;
   unit_cell.atom[2].z = 0;
   unit_cell.atom[2].mat = uc::internal::sublattice_materials ? 2 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[2].lc = 2;
   unit_cell.atom[2].hc = 0;
   unit_cell.atom[2].ni = 6;
   //---------------------------------------------
   unit_cell.atom[3].x = 0.25;
   unit_cell.atom[3].y = 0.25;
   unit_cell.atom[3].z = 0;
   unit_cell.atom[3].mat = uc::internal::sublattice_materials ? 3 : 0; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[3].lc = 0;
   unit_cell.atom[3].hc = 0;
   unit_cell.atom[3].ni = 6;
   //---------------------------------------------
   unit_cell.atom[4].x = 0;
   unit_cell.atom[4].y = 0;
   unit_cell.atom[4].z = 0.25;
   unit_cell.atom[4].mat = uc::internal::sublattice_materials ? 4 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[4].lc = 2;
   unit_cell.atom[4].hc = 1;
   unit_cell.atom[4].ni = 6;
   //---------------------------------------------
   unit_cell.atom[5].x = 0.25;
   unit_cell.atom[5].y = 0;
   unit_cell.atom[5].z = 0.25;
   unit_cell.atom[5].mat = uc::internal::sublattice_materials ? 5 : 0; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[5].lc = 0;
   unit_cell.atom[5].hc = 1;
   unit_cell.atom[5].ni = 6;
   //---------------------------------------------
   unit_cell.atom[6].x = 0;
   unit_cell.atom[6].y = 0.25;
   unit_cell.atom[6].z = 0.25;
   unit_cell.atom[6].mat = uc::internal::sublattice_materials ? 6 : 0; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[6].lc = 0;
   unit_cell.atom[6].hc = 1;
   unit_cell.atom[6].ni = 6;
   //---------------------------------------------
   unit_cell.atom[7].x = 0.25;
   unit_cell.atom[7].y = 0.25;
   unit_cell.atom[7].z = 0.25;
   unit_cell.atom[7].mat = uc::internal::sublattice_materials ? 7 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[7].lc = 2;
   unit_cell.atom[7].hc = 1;
   unit_cell.atom[7].ni = 6;
   //---------------------------------------------
   unit_cell.atom[8].x = 0.5;
   unit_cell.atom[8].y = 0;
   unit_cell.atom[8].z = 0;
   unit_cell.atom[8].mat = uc::internal::sublattice_materials ? 8 : 0; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[8].lc = 0;
   unit_cell.atom[8].hc = 0;
   unit_cell.atom[8].ni = 6;
   //---------------------------------------------
   unit_cell.atom[9].x = 0.75;
   unit_cell.atom[9].y = 0;
   unit_cell.atom[9].z = 0;
   unit_cell.atom[9].mat = uc::internal::sublattice_materials ? 9 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[9].lc = 2;
   unit_cell.atom[9].hc = 0;
   unit_cell.atom[9].ni = 6;
   //---------------------------------------------
   unit_cell.atom[10].x = 0.5;
   unit_cell.atom[10].y = 0.25;
   unit_cell.atom[10].z = 0;
   unit_cell.atom[10].mat = uc::internal::sublattice_materials ? 10 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[10].lc = 2;
   unit_cell.atom[10].hc = 0;
   unit_cell.atom[10].ni = 6;
   //---------------------------------------------
   unit_cell.atom[11].x = 0.75;
   unit_cell.atom[11].y = 0.25;
   unit_cell.atom[11].z = 0;
   unit_cell.atom[11].mat = uc::internal::sublattice_materials ? 11 : 1; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[11].lc = 1;
   unit_cell.atom[11].hc = 0;
   unit_cell.atom[11].ni = 6;
   //---------------------------------------------
   unit_cell.atom[12].x = 0.5;
   unit_cell.atom[12].y = 0;
   unit_cell.atom[12].z = 0.25;
   unit_cell.atom[12].mat = uc::internal::sublattice_materials ? 12 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[12].lc = 2;
   unit_cell.atom[12].hc = 1;
   unit_cell.atom[12].ni = 6;
   //---------------------------------------------
   unit_cell.atom[13].x = 0.75;
   unit_cell.atom[13].y = 0;
   unit_cell.atom[13].z = 0.25;
   unit_cell.atom[13].mat = uc::internal::sublattice_materials ? 13 : 1; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[13].lc = 1;
   unit_cell.atom[13].hc = 1;
   unit_cell.atom[13].ni = 6;
   //---------------------------------------------
   unit_cell.atom[14].x = 0.5;
   unit_cell.atom[14].y = 0.25;
   unit_cell.atom[14].z = 0.25;
   unit_cell.atom[14].mat = uc::internal::sublattice_materials ? 14 : 1; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[14].lc = 1;
   unit_cell.atom[14].hc = 1;
   unit_cell.atom[14].ni = 6;
   //---------------------------------------------
   unit_cell.atom[15].x = 0.75;
   unit_cell.atom[15].y = 0.25;
   unit_cell.atom[15].z = 0.25;
   unit_cell.atom[15].mat = uc::internal::sublattice_materials ? 15 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[15].lc = 2;
   unit_cell.atom[15].hc = 1;
   unit_cell.atom[15].ni = 6;
   //---------------------------------------------
   unit_cell.atom[16].x = 0;
   unit_cell.atom[16].y = 0.5;
   unit_cell.atom[16].z = 0;
   unit_cell.atom[16].mat = uc::internal::sublattice_materials ? 16 : 0; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[16].lc = 0;
   unit_cell.atom[16].hc = 0;
   unit_cell.atom[16].ni = 6;
   //---------------------------------------------
   unit_cell.atom[17].x = 0.25;
   unit_cell.atom[17].y = 0.5;
   unit_cell.atom[17].z = 0;
   unit_cell.atom[17].mat = uc::internal::sublattice_materials ? 17 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[17].lc = 2;
   unit_cell.atom[17].hc = 0;
   unit_cell.atom[17].ni = 6;
   //---------------------------------------------
   unit_cell.atom[18].x = 0;
   unit_cell.atom[18].y = 0.75;
   unit_cell.atom[18].z = 0;
   unit_cell.atom[18].mat = uc::internal::sublattice_materials ? 18 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[18].lc = 2;
   unit_cell.atom[18].hc = 0;
   unit_cell.atom[18].ni = 6;
   //---------------------------------------------
   unit_cell.atom[19].x = 0.25;
   unit_cell.atom[19].y = 0.75;
   unit_cell.atom[19].z = 0;
   unit_cell.atom[19].mat = uc::internal::sublattice_materials ? 19 : 1; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[19].lc = 1;
   unit_cell.atom[19].hc = 0;
   unit_cell.atom[19].ni = 6;
   //---------------------------------------------
   unit_cell.atom[20].x = 0;
   unit_cell.atom[20].y = 0.5;
   unit_cell.atom[20].z = 0.25;
   unit_cell.atom[20].mat = uc::internal::sublattice_materials ? 20 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[20].lc = 2;
   unit_cell.atom[20].hc = 1;
   unit_cell.atom[20].ni = 6;
   //---------------------------------------------
   unit_cell.atom[21].x = 0.25;
   unit_cell.atom[21].y = 0.5;
   unit_cell.atom[21].z = 0.25;
   unit_cell.atom[21].mat = uc::internal::sublattice_materials ? 21 : 1; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[21].lc = 1;
   unit_cell.atom[21].hc = 1;
   unit_cell.atom[21].ni = 6;
   //---------------------------------------------
   unit_cell.atom[22].x = 0;
   unit_cell.atom[22].y = 0.75;
   unit_cell.atom[22].z = 0.25;
   unit_cell.atom[22].mat = uc::internal::sublattice_materials ? 22 : 1; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[22].lc = 1;
   unit_cell.atom[22].hc = 1;
   unit_cell.atom[22].ni = 6;
   //---------------------------------------------
   unit_cell.atom[23].x = 0.25;
   unit_cell.atom[23].y = 0.75;
   unit_cell.atom[23].z = 0.25;
   unit_cell.atom[23].mat = uc::internal::sublattice_materials ? 23 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[23].lc = 2;
   unit_cell.atom[23].hc = 1;
   unit_cell.atom[23].ni = 6;
   //---------------------------------------------
   unit_cell.atom[24].x = 0.5;
   unit_cell.atom[24].y = 0.5;
   unit_cell.atom[24].z = 0;
   unit_cell.atom[24].mat = uc::internal::sublattice_materials ? 24 : 1; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[24].lc = 1;
   unit_cell.atom[24].hc = 0;
   unit_cell.atom[24].ni = 6;
   //---------------------------------------------
   unit_cell.atom[25].x = 0.75;
   unit_cell.atom[25].y = 0.5;
   unit_cell.atom[25].z = 0;
   unit_cell.atom[25].mat = uc::internal::sublattice_materials ? 25 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[25].lc = 2;
   unit_cell.atom[25].hc = 0;
   unit_cell.atom[25].ni = 6;
   //---------------------------------------------
   unit_cell.atom[26].x = 0.5;
   unit_cell.atom[26].y = 0.75;
   unit_cell.atom[26].z = 0;
   unit_cell.atom[26].mat = uc::internal::sublattice_materials ? 26 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[26].lc = 2;
   unit_cell.atom[26].hc = 0;
   unit_cell.atom[26].ni = 6;
   //---------------------------------------------
   unit_cell.atom[27].x = 0.75;
   unit_cell.atom[27].y = 0.75;
   unit_cell.atom[27].z = 0;
   unit_cell.atom[27].mat = uc::internal::sublattice_materials ? 27 : 0; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[27].lc = 0;
   unit_cell.atom[27].hc = 0;
   unit_cell.atom[27].ni = 6;
   //---------------------------------------------
   unit_cell.atom[28].x = 0.5;
   unit_cell.atom[28].y = 0.5;
   unit_cell.atom[28].z = 0.25;
   unit_cell.atom[28].mat = uc::internal::sublattice_materials ? 28 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[28].lc = 2;
   unit_cell.atom[28].hc = 1;
   unit_cell.atom[28].ni = 6;
   //---------------------------------------------
   unit_cell.atom[29].x = 0.75;
   unit_cell.atom[29].y = 0.5;
   unit_cell.atom[29].z = 0.25;
   unit_cell.atom[29].mat = uc::internal::sublattice_materials ? 29 : 0; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[29].lc = 0;
   unit_cell.atom[29].hc = 1;
   unit_cell.atom[29].ni = 6;
   //---------------------------------------------
   unit_cell.atom[30].x = 0.5;
   unit_cell.atom[30].y = 0.75;
   unit_cell.atom[30].z = 0.25;
   unit_cell.atom[30].mat = uc::internal::sublattice_materials ? 30 : 0; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[30].lc = 0;
   unit_cell.atom[30].hc = 1;
   unit_cell.atom[30].ni = 6;
   //---------------------------------------------
   unit_cell.atom[31].x = 0.75;
   unit_cell.atom[31].y = 0.75;
   unit_cell.atom[31].z = 0.25;
   unit_cell.atom[31].mat = uc::internal::sublattice_materials ? 31 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[31].lc = 2;
   unit_cell.atom[31].hc = 1;
   unit_cell.atom[31].ni = 6;
   //---------------------------------------------
   unit_cell.atom[32].x = 0;
   unit_cell.atom[32].y = 0;
   unit_cell.atom[32].z = 0.5;
   unit_cell.atom[32].mat = uc::internal::sublattice_materials ? 32 : 0; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[32].lc = 0;
   unit_cell.atom[32].hc = 2;
   unit_cell.atom[32].ni = 6;
   //---------------------------------------------
   unit_cell.atom[33].x = 0.25;
   unit_cell.atom[33].y = 0;
   unit_cell.atom[33].z = 0.5;
   unit_cell.atom[33].mat = uc::internal::sublattice_materials ? 33 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[33].lc = 2;
   unit_cell.atom[33].hc = 2;
   unit_cell.atom[33].ni = 6;
   //---------------------------------------------
   unit_cell.atom[34].x = 0;
   unit_cell.atom[34].y = 0.25;
   unit_cell.atom[34].z = 0.5;
   unit_cell.atom[34].mat = uc::internal::sublattice_materials ? 34 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[34].lc = 2;
   unit_cell.atom[34].hc = 2;
   unit_cell.atom[34].ni = 6;
   //---------------------------------------------
   unit_cell.atom[35].x = 0.25;
   unit_cell.atom[35].y = 0.25;
   unit_cell.atom[35].z = 0.5;
   unit_cell.atom[35].mat = uc::internal::sublattice_materials ? 35 : 1; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[35].lc = 1;
   unit_cell.atom[35].hc = 2;
   unit_cell.atom[35].ni = 6;
   //---------------------------------------------
   unit_cell.atom[36].x = 0;
   unit_cell.atom[36].y = 0;
   unit_cell.atom[36].z = 0.75;
   unit_cell.atom[36].mat = uc::internal::sublattice_materials ? 36 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[36].lc = 2;
   unit_cell.atom[36].hc = 3;
   unit_cell.atom[36].ni = 6;
   //---------------------------------------------
   unit_cell.atom[37].x = 0.25;
   unit_cell.atom[37].y = 0;
   unit_cell.atom[37].z = 0.75;
   unit_cell.atom[37].mat = uc::internal::sublattice_materials ? 37 : 1; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[37].lc = 1;
   unit_cell.atom[37].hc = 3;
   unit_cell.atom[37].ni = 6;
   //---------------------------------------------
   unit_cell.atom[38].x = 0;
   unit_cell.atom[38].y = 0.25;
   unit_cell.atom[38].z = 0.75;
   unit_cell.atom[38].mat = uc::internal::sublattice_materials ? 38 : 1; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[38].lc = 1;
   unit_cell.atom[38].hc = 3;
   unit_cell.atom[38].ni = 6;
   //---------------------------------------------
   unit_cell.atom[39].x = 0.25;
   unit_cell.atom[39].y = 0.25;
   unit_cell.atom[39].z = 0.75;
   unit_cell.atom[39].mat = uc::internal::sublattice_materials ? 39 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[39].lc = 2;
   unit_cell.atom[39].hc = 3;
   unit_cell.atom[39].ni = 6;
   //---------------------------------------------
   unit_cell.atom[40].x = 0.5;
   unit_cell.atom[40].y = 0;
   unit_cell.atom[40].z = 0.5;
   unit_cell.atom[40].mat = uc::internal::sublattice_materials ? 40 : 1; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[40].lc = 1;
   unit_cell.atom[40].hc = 2;
   unit_cell.atom[40].ni = 6;
   //---------------------------------------------
   unit_cell.atom[41].x = 0.75;
   unit_cell.atom[41].y = 0;
   unit_cell.atom[41].z = 0.5;
   unit_cell.atom[41].mat = uc::internal::sublattice_materials ? 41 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[41].lc = 2;
   unit_cell.atom[41].hc = 2;
   unit_cell.atom[41].ni = 6;
   //---------------------------------------------
   unit_cell.atom[42].x = 0.5;
   unit_cell.atom[42].y = 0.25;
   unit_cell.atom[42].z = 0.5;
   unit_cell.atom[42].mat = uc::internal::sublattice_materials ? 42 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[42].lc = 2;
   unit_cell.atom[42].hc = 2;
   unit_cell.atom[42].ni = 6;
   //---------------------------------------------
   unit_cell.atom[43].x = 0.75;
   unit_cell.atom[43].y = 0.25;
   unit_cell.atom[43].z = 0.5;
   unit_cell.atom[43].mat = uc::internal::sublattice_materials ? 43 : 0; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[43].lc = 0;
   unit_cell.atom[43].hc = 2;
   unit_cell.atom[43].ni = 6;
   //---------------------------------------------
   unit_cell.atom[44].x = 0.5;
   unit_cell.atom[44].y = 0;
   unit_cell.atom[44].z = 0.75;
   unit_cell.atom[44].mat = uc::internal::sublattice_materials ? 44 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[44].lc = 2;
   unit_cell.atom[44].hc = 3;
   unit_cell.atom[44].ni = 6;
   //---------------------------------------------
   unit_cell.atom[45].x = 0.75;
   unit_cell.atom[45].y = 0;
   unit_cell.atom[45].z = 0.75;
   unit_cell.atom[45].mat = uc::internal::sublattice_materials ? 45 : 0; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[45].lc = 0;
   unit_cell.atom[45].hc = 3;
   unit_cell.atom[45].ni = 6;
   //---------------------------------------------
   unit_cell.atom[46].x = 0.5;
   unit_cell.atom[46].y = 0.25;
   unit_cell.atom[46].z = 0.75;
   unit_cell.atom[46].mat = uc::internal::sublattice_materials ? 46 : 0; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[46].lc = 0;
   unit_cell.atom[46].hc = 3;
   unit_cell.atom[46].ni = 6;
   //---------------------------------------------
   unit_cell.atom[47].x = 0.75;
   unit_cell.atom[47].y = 0.25;
   unit_cell.atom[47].z = 0.75;
   unit_cell.atom[47].mat = uc::internal::sublattice_materials ? 47 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[47].lc = 2;
   unit_cell.atom[47].hc = 3;
   unit_cell.atom[47].ni = 6;
   //---------------------------------------------
   unit_cell.atom[48].x = 0;
   unit_cell.atom[48].y = 0.5;
   unit_cell.atom[48].z = 0.5;
   unit_cell.atom[48].mat = uc::internal::sublattice_materials ? 48 : 1; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[48].lc = 1;
   unit_cell.atom[48].hc = 2;
   unit_cell.atom[48].ni = 6;
   //---------------------------------------------
   unit_cell.atom[49].x = 0.25;
   unit_cell.atom[49].y = 0.5;
   unit_cell.atom[49].z = 0.5;
   unit_cell.atom[49].mat = uc::internal::sublattice_materials ? 49 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[49].lc = 2;
   unit_cell.atom[49].hc = 2;
   unit_cell.atom[49].ni = 6;
   //---------------------------------------------
   unit_cell.atom[50].x = 0;
   unit_cell.atom[50].y = 0.75;
   unit_cell.atom[50].z = 0.5;
   unit_cell.atom[50].mat = uc::internal::sublattice_materials ? 50 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[50].lc = 2;
   unit_cell.atom[50].hc = 2;
   unit_cell.atom[50].ni = 6;
   //---------------------------------------------
   unit_cell.atom[51].x = 0.25;
   unit_cell.atom[51].y = 0.75;
   unit_cell.atom[51].z = 0.5;
   unit_cell.atom[51].mat = uc::internal::sublattice_materials ? 51 : 0; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[51].lc = 0;
   unit_cell.atom[51].hc = 2;
   unit_cell.atom[51].ni = 6;
   //---------------------------------------------
   unit_cell.atom[52].x = 0;
   unit_cell.atom[52].y = 0.5;
   unit_cell.atom[52].z = 0.75;
   unit_cell.atom[52].mat = uc::internal::sublattice_materials ? 52 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[52].lc = 2;
   unit_cell.atom[52].hc = 3;
   unit_cell.atom[52].ni = 6;
   //---------------------------------------------
   unit_cell.atom[53].x = 0.25;
   unit_cell.atom[53].y = 0.5;
   unit_cell.atom[53].z = 0.75;
   unit_cell.atom[53].mat = uc::internal::sublattice_materials ? 53 : 0; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[53].lc = 0;
   unit_cell.atom[53].hc = 3;
   unit_cell.atom[53].ni = 6;
   //---------------------------------------------
   unit_cell.atom[54].x = 0;
   unit_cell.atom[54].y = 0.75;
   unit_cell.atom[54].z = 0.75;
   unit_cell.atom[54].mat = uc::internal::sublattice_materials ? 54 : 0; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[54].lc = 0;
   unit_cell.atom[54].hc = 3;
   unit_cell.atom[54].ni = 6;
   //---------------------------------------------
   unit_cell.atom[55].x = 0.25;
   unit_cell.atom[55].y = 0.75;
   unit_cell.atom[55].z = 0.75;
   unit_cell.atom[55].mat = uc::internal::sublattice_materials ? 55 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[55].lc = 2;
   unit_cell.atom[55].hc = 3;
   unit_cell.atom[55].ni = 6;
   //---------------------------------------------
   unit_cell.atom[56].x = 0.5;
   unit_cell.atom[56].y = 0.5;
   unit_cell.atom[56].z = 0.5;
   unit_cell.atom[56].mat = uc::internal::sublattice_materials ? 56 : 0; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[56].lc = 0;
   unit_cell.atom[56].hc = 2;
   unit_cell.atom[56].ni = 6;
   //---------------------------------------------
   unit_cell.atom[57].x = 0.75;
   unit_cell.atom[57].y = 0.5;
   unit_cell.atom[57].z = 0.5;
   unit_cell.atom[57].mat = uc::internal::sublattice_materials ? 57 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[57].lc = 2;
   unit_cell.atom[57].hc = 2;
   unit_cell.atom[57].ni = 6;
   //---------------------------------------------
   unit_cell.atom[58].x = 0.5;
   unit_cell.atom[58].y = 0.75;
   unit_cell.atom[58].z = 0.5;
   unit_cell.atom[58].mat = uc::internal::sublattice_materials ? 58 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[58].lc = 2;
   unit_cell.atom[58].hc = 2;
   unit_cell.atom[58].ni = 6;
   //---------------------------------------------
   unit_cell.atom[59].x = 0.75;
   unit_cell.atom[59].y = 0.75;
   unit_cell.atom[59].z = 0.5;
   unit_cell.atom[59].mat = uc::internal::sublattice_materials ? 59 : 1; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[59].lc = 1;
   unit_cell.atom[59].hc = 2;
   unit_cell.atom[59].ni = 6;
   //---------------------------------------------
   unit_cell.atom[60].x = 0.5;
   unit_cell.atom[60].y = 0.5;
   unit_cell.atom[60].z = 0.75;
   unit_cell.atom[60].mat = uc::internal::sublattice_materials ? 60 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[60].lc = 2;
   unit_cell.atom[60].hc = 3;
   unit_cell.atom[60].ni = 6;
   //---------------------------------------------
   unit_cell.atom[61].x = 0.75;
   unit_cell.atom[61].y = 0.5;
   unit_cell.atom[61].z = 0.75;
   unit_cell.atom[61].mat = uc::internal::sublattice_materials ? 61 : 1; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[61].lc = 1;
   unit_cell.atom[61].hc = 3;
   unit_cell.atom[61].ni = 6;
   //---------------------------------------------
   unit_cell.atom[62].x = 0.5;
   unit_cell.atom[62].y = 0.75;
   unit_cell.atom[62].z = 0.75;
   unit_cell.atom[62].mat = uc::internal::sublattice_materials ? 62 : 1; // if sublattice material is defined, then identify as metal sublattice
   unit_cell.atom[62].lc = 1;
   unit_cell.atom[62].hc = 3;
   unit_cell.atom[62].ni = 6;
   //---------------------------------------------
   unit_cell.atom[63].x = 0.75;
   unit_cell.atom[63].y = 0.75;
   unit_cell.atom[63].z = 0.75;
   unit_cell.atom[63].mat = uc::internal::sublattice_materials ? 63 : 2; // if sublattice material is defined, then identify as oxygen sublattice
   unit_cell.atom[63].lc = 2;
   unit_cell.atom[63].hc = 3;
   unit_cell.atom[63].ni = 6;
   //---------------------------------------------
   unit_cell.cutoff_radius = 0.51; // normalised to unit cell size (change to 0.26 for oxygen NNs)

   uc::internal::calculate_interactions(unit_cell);

   // Set actual unit cell size after calculating interactions
   unit_cell.dimensions[0] *= unitcell::internal::unit_cell_size_x;
   unit_cell.dimensions[1] *= unitcell::internal::unit_cell_size_y;
   unit_cell.dimensions[2] *= unitcell::internal::unit_cell_size_z;

   std::cout << "Generated rocksalt-supercell with " << unit_cell.atom.size() << " atoms" << std::endl;

   return;

}

} // end of internal namespace
} // end of unitcell namespace
