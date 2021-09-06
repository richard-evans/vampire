// This program can be used to output data for plotting radial distribution functions.
// This program is run on its own, i.e. not compiled with vampire.
#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include<iomanip>

// Make sure to use correct file path to include unitcell.hpp
#include"/Users/jbc525/Documents/vampire/hdr/unitcell.hpp"

// Note: includes code modified from interactions.cpp

//------------------------------------------------------------------------------
// Class and struct definitions
//------------------------------------------------------------------------------

namespace local{

   class atom_t{

   public:
      int mat; // material
      int lc;
      int hc;
      double x; // positions
      double y;
      double z;
      int id; // atom number in unit cell
      int idx; // unit cell number
      int idy;
      int idz;
      bool nm;
   };

// Class to hold information for interactions between central unit cell atoms 
// and all other atoms. Values held in instances of this class can be used to 
// plot radial distribution functions
   struct interaction_t{
      unsigned int mat; /// Material
      unsigned int lc; /// Lattice category
      unsigned int hc; /// Height category
      unsigned int inter_mat; /// Interactor material
      unsigned int inter_hc; /// Interactor height category
      unsigned int bin; /// Radial bin interactor resides in

      //constructor
      interaction_t():
         mat(0),
         lc(0),
         hc(0),
         inter_mat(0),
         inter_hc(0),
         bin(0)
      {
      };
   };

} // End local namespace

//------------------------------------------------------------------------------
// Function definitions
//------------------------------------------------------------------------------

// This function initialises a unit cell instance of NdFeB, you can replace this with your desired unitcell.

void build_NdFeB_Cell(unitcell::unit_cell_t& unit_cell){

   // Set basic unit cell properties
   unit_cell.dimensions[0] = 1.0;
   unit_cell.dimensions[1] = 1.0;
   unit_cell.dimensions[2] = 12.2 / 8.8;

   unit_cell.shape[0][0] = 1;
   unit_cell.shape[0][1] = 0;
   unit_cell.shape[0][2] = 0;

   unit_cell.shape[1][0] = 0;
   unit_cell.shape[1][1] = 1;
   unit_cell.shape[1][2] = 0;

   unit_cell.shape[2][0] = 0;
   unit_cell.shape[2][1] = 0;
   unit_cell.shape[2][2] = 1;

   unit_cell.lcsize = 68;
   unit_cell.hcsize = 4;
   unit_cell.interaction_range = 1;
   unit_cell.atom.resize(68);

   //-----------------------------
   unit_cell.atom[0].x   = 0.268;
   unit_cell.atom[0].y   = 0.268;
   unit_cell.atom[0].z   = 0;
   unit_cell.atom[0].mat = 0;
   unit_cell.atom[0].lc  = 0;
   unit_cell.atom[0].hc  = 0;
   unit_cell.atom[0].nm  = false;
   //-----------------------------
   unit_cell.atom[1].x   = 0.732;
   unit_cell.atom[1].y   = 0.732;
   unit_cell.atom[1].z   = 0;
   unit_cell.atom[1].mat = 0;
   unit_cell.atom[1].lc  = 1;
   unit_cell.atom[1].hc  = 0;
   unit_cell.atom[1].nm  = false;
   //-----------------------------
   unit_cell.atom[2].x   = 0.232;
   unit_cell.atom[2].y   = 0.768;
   unit_cell.atom[2].z   = 0.5;
   unit_cell.atom[2].mat = 0;
   unit_cell.atom[2].lc  = 2;
   unit_cell.atom[2].hc  = 2;
   unit_cell.atom[2].nm  = false;
   //-----------------------------
   unit_cell.atom[3].x   = 0.768;
   unit_cell.atom[3].y   = 0.232;
   unit_cell.atom[3].z   = 0.5;
   unit_cell.atom[3].mat = 0;
   unit_cell.atom[3].lc  = 3;
   unit_cell.atom[3].hc  = 2;
   unit_cell.atom[3].nm  = false;
   //-----------------------------
   unit_cell.atom[4].x   = 0.14;
   unit_cell.atom[4].y   = 0.86;
   unit_cell.atom[4].z   = 0;
   unit_cell.atom[4].mat = 0;
   unit_cell.atom[4].lc  = 4;
   unit_cell.atom[4].hc  = 0;
   unit_cell.atom[4].nm  = false;
   //-----------------------------
   unit_cell.atom[5].x   = 0.86;
   unit_cell.atom[5].y   = 0.14;
   unit_cell.atom[5].z   = 0;
   unit_cell.atom[5].mat = 0;
   unit_cell.atom[5].lc  = 5;
   unit_cell.atom[5].hc  = 0;
   unit_cell.atom[5].nm  = false;
   //-----------------------------
   unit_cell.atom[6].x   = 0.64;
   unit_cell.atom[6].y   = 0.64;
   unit_cell.atom[6].z   = 0.5;
   unit_cell.atom[6].mat = 0;
   unit_cell.atom[6].lc  = 6;
   unit_cell.atom[6].hc  = 2;
   unit_cell.atom[6].nm  = false;
   //-----------------------------
   unit_cell.atom[7].x   = 0.36;
   unit_cell.atom[7].y   = 0.36;
   unit_cell.atom[7].z   = 0.5;
   unit_cell.atom[7].mat = 0;
   unit_cell.atom[7].lc  = 7;
   unit_cell.atom[7].hc  = 2;
   unit_cell.atom[7].nm  = false;
   //-----------------------------
   unit_cell.atom[8].x   = 0.223;
   unit_cell.atom[8].y   = 0.567;
   unit_cell.atom[8].z   = 0.127;
   unit_cell.atom[8].mat = 1;
   unit_cell.atom[8].lc  = 8;
   unit_cell.atom[8].hc  = 1;
   unit_cell.atom[8].nm  = false;
   //-----------------------------
   unit_cell.atom[9].x   = 0.777;
   unit_cell.atom[9].y   = 0.433;
   unit_cell.atom[9].z   = 0.127;
   unit_cell.atom[9].mat = 1;
   unit_cell.atom[9].lc  = 9;
   unit_cell.atom[9].hc  = 1;
   unit_cell.atom[9].nm  = false;
   //-----------------------------
   unit_cell.atom[10].x   = 0.933;
   unit_cell.atom[10].y   = 0.723;
   unit_cell.atom[10].z   = 0.627;
   unit_cell.atom[10].mat = 1;
   unit_cell.atom[10].lc  = 10;
   unit_cell.atom[10].hc  = 3;
   unit_cell.atom[10].nm  = false;
   //-----------------------------
   unit_cell.atom[11].x   = 0.067;
   unit_cell.atom[11].y   = 0.277;
   unit_cell.atom[11].z   = 0.627;
   unit_cell.atom[11].mat = 1;
   unit_cell.atom[11].lc  = 11;
   unit_cell.atom[11].hc  = 3;
   unit_cell.atom[11].nm  = false;
   //-----------------------------
   unit_cell.atom[12].x   = 0.277;
   unit_cell.atom[12].y   = 0.067;
   unit_cell.atom[12].z   = 0.373;
   unit_cell.atom[12].mat = 1;
   unit_cell.atom[12].lc  = 12;
   unit_cell.atom[12].hc  = 1;
   unit_cell.atom[12].nm  = false;
   //-----------------------------
   unit_cell.atom[13].x   = 0.723;
   unit_cell.atom[13].y   = 0.933;
   unit_cell.atom[13].z   = 0.373;
   unit_cell.atom[13].mat = 1;
   unit_cell.atom[13].lc  = 13;
   unit_cell.atom[13].hc  = 1;
   unit_cell.atom[13].nm  = false;
   //-----------------------------
   unit_cell.atom[14].x   = 0.567;
   unit_cell.atom[14].y   = 0.223;
   unit_cell.atom[14].z   = 0.873;
   unit_cell.atom[14].mat = 1;
   unit_cell.atom[14].lc  = 14;
   unit_cell.atom[14].hc  = 3;
   unit_cell.atom[14].nm  = false;
   //-----------------------------
   unit_cell.atom[15].x   = 0.433;
   unit_cell.atom[15].y   = 0.777;
   unit_cell.atom[15].z   = 0.873;
   unit_cell.atom[15].mat = 1;
   unit_cell.atom[15].lc  = 15;
   unit_cell.atom[15].hc  = 3;
   unit_cell.atom[15].nm  = false;
   //-----------------------------
   unit_cell.atom[16].x   = 0.777;
   unit_cell.atom[16].y   = 0.433;
   unit_cell.atom[16].z   = 0.873;
   unit_cell.atom[16].mat = 1;
   unit_cell.atom[16].lc  = 16;
   unit_cell.atom[16].hc  = 3;
   unit_cell.atom[16].nm  = false;
   //-----------------------------
   unit_cell.atom[17].x   = 0.223;
   unit_cell.atom[17].y   = 0.567;
   unit_cell.atom[17].z   = 0.873;
   unit_cell.atom[17].mat = 1;
   unit_cell.atom[17].lc  = 17;
   unit_cell.atom[17].hc  = 3;
   unit_cell.atom[17].nm  = false; 
   //-----------------------------
   unit_cell.atom[18].x   = 0.067;
   unit_cell.atom[18].y   = 0.277;
   unit_cell.atom[18].z   = 0.373;
   unit_cell.atom[18].mat = 1;
   unit_cell.atom[18].lc  = 18;
   unit_cell.atom[18].hc  = 1;
   unit_cell.atom[18].nm  = false;
   //-----------------------------
   unit_cell.atom[19].x   = 0.933;
   unit_cell.atom[19].y   = 0.723;
   unit_cell.atom[19].z   = 0.373;
   unit_cell.atom[19].mat = 1;
   unit_cell.atom[19].lc  = 19;
   unit_cell.atom[19].hc  = 1;
   unit_cell.atom[19].nm  = false;
   //-----------------------------
   unit_cell.atom[20].x   = 0.723;
   unit_cell.atom[20].y   = 0.933;
   unit_cell.atom[20].z   = 0.627;
   unit_cell.atom[20].mat = 1;
   unit_cell.atom[20].lc  = 20;
   unit_cell.atom[20].hc  = 3;
   unit_cell.atom[20].nm  = false;
   //-----------------------------
   unit_cell.atom[21].x   = 0.277;
   unit_cell.atom[21].y   = 0.067;
   unit_cell.atom[21].z   = 0.627;
   unit_cell.atom[21].mat = 1;
   unit_cell.atom[21].lc  = 21;
   unit_cell.atom[21].hc  = 3;
   unit_cell.atom[21].nm  = false;
   //-----------------------------
   unit_cell.atom[22].x   = 0.433;
   unit_cell.atom[22].y   = 0.777;
   unit_cell.atom[22].z   = 0.127;
   unit_cell.atom[22].mat = 1;
   unit_cell.atom[22].lc  = 22;
   unit_cell.atom[22].hc  = 1;
   unit_cell.atom[22].nm  = false;
   //-----------------------------
   unit_cell.atom[23].x   = 0.567;
   unit_cell.atom[23].y   = 0.223;
   unit_cell.atom[23].z   = 0.127;
   unit_cell.atom[23].mat = 1;
   unit_cell.atom[23].lc  = 23;
   unit_cell.atom[23].hc  = 1;
   unit_cell.atom[23].nm  = false;
   //-----------------------------
   unit_cell.atom[24].x   = 0.037;
   unit_cell.atom[24].y   = 0.36;
   unit_cell.atom[24].z   = 0.176;
   unit_cell.atom[24].mat = 1;
   unit_cell.atom[24].lc  = 24;
   unit_cell.atom[24].hc  = 1;
   unit_cell.atom[24].nm  = false;
   //-----------------------------
   unit_cell.atom[25].x   = 0.963;
   unit_cell.atom[25].y   = 0.64;
   unit_cell.atom[25].z   = 0.176;
   unit_cell.atom[25].mat = 1;
   unit_cell.atom[25].lc  = 25;
   unit_cell.atom[25].hc  = 1;
   unit_cell.atom[25].nm  = false;
   //-----------------------------
   unit_cell.atom[26].x   = 0.14;
   unit_cell.atom[26].y   = 0.537;
   unit_cell.atom[26].z   = 0.676;
   unit_cell.atom[26].mat = 1;
   unit_cell.atom[26].lc  = 26;
   unit_cell.atom[26].hc  = 3;
   unit_cell.atom[26].nm  = false;
   //-----------------------------
   unit_cell.atom[27].x   = 0.86;
   unit_cell.atom[27].y   = 0.463;
   unit_cell.atom[27].z   = 0.676;
   unit_cell.atom[27].mat = 1;
   unit_cell.atom[27].lc  = 27;
   unit_cell.atom[27].hc  = 3;
   unit_cell.atom[27].nm  = false;
   //-----------------------------
   unit_cell.atom[28].x   = 0.463;
   unit_cell.atom[28].y   = 0.86;
   unit_cell.atom[28].z   = 0.324;
   unit_cell.atom[28].mat = 1;
   unit_cell.atom[28].lc  = 28;
   unit_cell.atom[28].hc  = 1;
   unit_cell.atom[28].nm  = false;
   //-----------------------------
   unit_cell.atom[29].x   = 0.537;
   unit_cell.atom[29].y   = 0.14;
   unit_cell.atom[29].z   = 0.324;
   unit_cell.atom[29].mat = 1;
   unit_cell.atom[29].lc  = 29;
   unit_cell.atom[29].hc  = 1;
   unit_cell.atom[29].nm  = false;
   //-----------------------------
   unit_cell.atom[30].x   = 0.36;
   unit_cell.atom[30].y   = 0.037;
   unit_cell.atom[30].z   = 0.824;
   unit_cell.atom[30].mat = 1;
   unit_cell.atom[30].lc  = 30;
   unit_cell.atom[30].hc  = 3;
   unit_cell.atom[30].nm  = false;
   //-----------------------------
   unit_cell.atom[31].x   = 0.64;
   unit_cell.atom[31].y   = 0.963;
   unit_cell.atom[31].z   = 0.824;
   unit_cell.atom[31].mat = 1;
   unit_cell.atom[31].lc  = 31;
   unit_cell.atom[31].hc  = 3;
   unit_cell.atom[31].nm  = false;
   //-----------------------------
   unit_cell.atom[32].x   = 0.963;
   unit_cell.atom[32].y   = 0.64;
   unit_cell.atom[32].z   = 0.824;
   unit_cell.atom[32].mat = 1;
   unit_cell.atom[32].lc  = 32;
   unit_cell.atom[32].hc  = 3;
   unit_cell.atom[32].nm  = false;   
   //-----------------------------
   unit_cell.atom[33].x   = 0.037;
   unit_cell.atom[33].y   = 0.36;
   unit_cell.atom[33].z   = 0.824;
   unit_cell.atom[33].mat = 1;
   unit_cell.atom[33].lc  = 33;
   unit_cell.atom[33].hc  = 3;
   unit_cell.atom[33].nm  = false;
   //-----------------------------
   unit_cell.atom[34].x   = 0.86;
   unit_cell.atom[34].y   = 0.463;
   unit_cell.atom[34].z   = 0.324;
   unit_cell.atom[34].mat = 1;
   unit_cell.atom[34].lc  = 34;
   unit_cell.atom[34].hc  = 1;
   unit_cell.atom[34].nm  = false;
   //-----------------------------
   unit_cell.atom[35].x   = 0.14;
   unit_cell.atom[35].y   = 0.537;
   unit_cell.atom[35].z   = 0.324;
   unit_cell.atom[35].mat = 1;
   unit_cell.atom[35].lc  = 35;
   unit_cell.atom[35].hc  = 1;
   unit_cell.atom[35].nm  = false;
   //-----------------------------
   unit_cell.atom[36].x   = 0.537;
   unit_cell.atom[36].y   = 0.14;
   unit_cell.atom[36].z   = 0.676;
   unit_cell.atom[36].mat = 1;
   unit_cell.atom[36].lc  = 36;
   unit_cell.atom[36].hc  = 3;
   unit_cell.atom[36].nm  = false;
   //-----------------------------
   unit_cell.atom[37].x   = 0.463;
   unit_cell.atom[37].y   = 0.86;
   unit_cell.atom[37].z   = 0.676;
   unit_cell.atom[37].mat = 1;
   unit_cell.atom[37].lc  = 37;
   unit_cell.atom[37].hc  = 3;
   unit_cell.atom[37].nm  = false;
   //-----------------------------
   unit_cell.atom[38].x   = 0.64;
   unit_cell.atom[38].y   = 0.963;
   unit_cell.atom[38].z   = 0.176;
   unit_cell.atom[38].mat = 1;
   unit_cell.atom[38].lc  = 38;
   unit_cell.atom[38].hc  = 1;
   unit_cell.atom[38].nm  = false;
   //-----------------------------
   unit_cell.atom[39].x   = 0.36;
   unit_cell.atom[39].y   = 0.037;
   unit_cell.atom[39].z   = 0.176;
   unit_cell.atom[39].mat = 1;
   unit_cell.atom[39].lc  = 39;
   unit_cell.atom[39].hc  = 1;
   unit_cell.atom[39].nm  = false;
   //-----------------------------
   unit_cell.atom[40].x   = 0.098;
   unit_cell.atom[40].y   = 0.098;
   unit_cell.atom[40].z   = 0.204;
   unit_cell.atom[40].mat = 1;
   unit_cell.atom[40].lc  = 40;
   unit_cell.atom[40].hc  = 1;
   unit_cell.atom[40].nm  = false;
   //-----------------------------
   unit_cell.atom[41].x   = 0.902;
   unit_cell.atom[41].y   = 0.902;
   unit_cell.atom[41].z   = 0.204;
   unit_cell.atom[41].mat = 1;
   unit_cell.atom[41].lc  = 41;
   unit_cell.atom[41].hc  = 1;
   unit_cell.atom[41].nm  = false;
   //-----------------------------
   unit_cell.atom[42].x   = 0.402;
   unit_cell.atom[42].y   = 0.598;
   unit_cell.atom[42].z   = 0.704;
   unit_cell.atom[42].mat = 1;
   unit_cell.atom[42].lc  = 42;
   unit_cell.atom[42].hc  = 3;
   unit_cell.atom[42].nm  = false;
   //-----------------------------
   unit_cell.atom[43].x   = 0.598;
   unit_cell.atom[43].y   = 0.402;
   unit_cell.atom[43].z   = 0.704;
   unit_cell.atom[43].mat = 1;
   unit_cell.atom[43].lc  = 43;
   unit_cell.atom[43].hc  = 3;
   unit_cell.atom[43].nm  = false;
   //-----------------------------
   unit_cell.atom[44].x   = 0.402;
   unit_cell.atom[44].y   = 0.598;
   unit_cell.atom[44].z   = 0.296;
   unit_cell.atom[44].mat = 1;
   unit_cell.atom[44].lc  = 44;
   unit_cell.atom[44].hc  = 1;
   unit_cell.atom[44].nm  = false;
   //-----------------------------
   unit_cell.atom[45].x   = 0.598;
   unit_cell.atom[45].y   = 0.402;
   unit_cell.atom[45].z   = 0.296;
   unit_cell.atom[45].mat = 1;
   unit_cell.atom[45].lc  = 45;
   unit_cell.atom[45].hc  = 1;
   unit_cell.atom[45].nm  = false;
   //-----------------------------
   unit_cell.atom[46].x   = 0.098;
   unit_cell.atom[46].y   = 0.098;
   unit_cell.atom[46].z   = 0.796;
   unit_cell.atom[46].mat = 1;
   unit_cell.atom[46].lc  = 46;
   unit_cell.atom[46].hc  = 3;
   unit_cell.atom[46].nm  = false;
   //-----------------------------
   unit_cell.atom[47].x   = 0.902;
   unit_cell.atom[47].y   = 0.902;
   unit_cell.atom[47].z   = 0.796;
   unit_cell.atom[47].mat = 1;
   unit_cell.atom[47].lc  = 47;
   unit_cell.atom[47].hc  = 3;
   unit_cell.atom[47].nm  = false;
   //-----------------------------
   unit_cell.atom[48].x   = 0.317;
   unit_cell.atom[48].y   = 0.317;
   unit_cell.atom[48].z   = 0.246;
   unit_cell.atom[48].mat = 1;
   unit_cell.atom[48].lc  = 48;
   unit_cell.atom[48].hc  = 1;
   unit_cell.atom[48].nm  = false;
   //-----------------------------
   unit_cell.atom[49].x   = 0.683;
   unit_cell.atom[49].y   = 0.683;
   unit_cell.atom[49].z   = 0.246;
   unit_cell.atom[49].mat = 1;
   unit_cell.atom[49].lc  = 49;
   unit_cell.atom[49].hc  = 1;
   unit_cell.atom[49].nm  = false;
   //-----------------------------
   unit_cell.atom[50].x   = 0.183;
   unit_cell.atom[50].y   = 0.817;
   unit_cell.atom[50].z   = 0.746;
   unit_cell.atom[50].mat = 1;
   unit_cell.atom[50].lc  = 50;
   unit_cell.atom[50].hc  = 3;
   unit_cell.atom[50].nm  = false;
   //-----------------------------
   unit_cell.atom[51].x   = 0.817;
   unit_cell.atom[51].y   = 0.183;
   unit_cell.atom[51].z   = 0.746;
   unit_cell.atom[51].mat = 1;
   unit_cell.atom[51].lc  = 51;
   unit_cell.atom[51].hc  = 3;
   unit_cell.atom[51].nm  = false;
   //-----------------------------
   unit_cell.atom[52].x   = 0.183;
   unit_cell.atom[52].y   = 0.817;
   unit_cell.atom[52].z   = 0.254;
   unit_cell.atom[52].mat = 1;
   unit_cell.atom[52].lc  = 52;
   unit_cell.atom[52].hc  = 1;
   unit_cell.atom[52].nm  = false;
   //-----------------------------
   unit_cell.atom[53].x   = 0.817;
   unit_cell.atom[53].y   = 0.183;
   unit_cell.atom[53].z   = 0.254;
   unit_cell.atom[53].mat = 1;
   unit_cell.atom[53].lc  = 53;
   unit_cell.atom[53].hc  = 1;
   unit_cell.atom[53].nm  = false;
   //-----------------------------
   unit_cell.atom[54].x   = 0.317;
   unit_cell.atom[54].y   = 0.317;
   unit_cell.atom[54].z   = 0.754;
   unit_cell.atom[54].mat = 1;
   unit_cell.atom[54].lc  = 54;
   unit_cell.atom[54].hc  = 3;
   unit_cell.atom[54].nm  = false;
   //-----------------------------
   unit_cell.atom[55].x   = 0.683;
   unit_cell.atom[55].y   = 0.683;
   unit_cell.atom[55].z   = 0.754;
   unit_cell.atom[55].mat = 1;
   unit_cell.atom[55].lc  = 55;
   unit_cell.atom[55].hc  = 3;
   unit_cell.atom[55].nm  = false;
   //-----------------------------
   unit_cell.atom[56].x   = 0.5;
   unit_cell.atom[56].y   = 0.5;
   unit_cell.atom[56].z   = 0.114;
   unit_cell.atom[56].mat = 1;
   unit_cell.atom[56].lc  = 56;
   unit_cell.atom[56].hc  = 1;
   unit_cell.atom[56].nm  = false;
   //-----------------------------
   unit_cell.atom[57].x   = 0;
   unit_cell.atom[57].y   = 0;
   unit_cell.atom[57].z   = 0.614;
   unit_cell.atom[57].mat = 1;
   unit_cell.atom[57].lc  = 57;
   unit_cell.atom[57].hc  = 3;
   unit_cell.atom[57].nm  = false;
   //-----------------------------
   unit_cell.atom[58].x   = 0;
   unit_cell.atom[58].y   = 0;
   unit_cell.atom[58].z   = 0.386;
   unit_cell.atom[58].mat = 1;
   unit_cell.atom[58].lc  = 58;
   unit_cell.atom[58].hc  = 1;
   unit_cell.atom[58].nm  = false;
   //-----------------------------
   unit_cell.atom[59].x   = 0.5;
   unit_cell.atom[59].y   = 0.5;
   unit_cell.atom[59].z   = 0.886;
   unit_cell.atom[59].mat = 1;
   unit_cell.atom[59].lc  = 59;
   unit_cell.atom[59].hc  = 3;
   unit_cell.atom[59].nm  = false;
   //-----------------------------
   unit_cell.atom[60].x   = 0;
   unit_cell.atom[60].y   = 0.5;
   unit_cell.atom[60].z   = 0;
   unit_cell.atom[60].mat = 1;
   unit_cell.atom[60].lc  = 60;
   unit_cell.atom[60].hc  = 0;
   unit_cell.atom[60].nm  = false;
   //-----------------------------
   unit_cell.atom[61].x   = 0;
   unit_cell.atom[61].y   = 0.5;
   unit_cell.atom[61].z   = 0.5;
   unit_cell.atom[61].mat = 1;
   unit_cell.atom[61].lc  = 61;
   unit_cell.atom[61].hc  = 2;
   unit_cell.atom[61].nm  = false;
   //-----------------------------
   unit_cell.atom[62].x   = 0.5;
   unit_cell.atom[62].y   = 0;
   unit_cell.atom[62].z   = 0.5;
   unit_cell.atom[62].mat = 1;
   unit_cell.atom[62].lc  = 62;
   unit_cell.atom[62].hc  = 2;
   unit_cell.atom[62].nm  = false;
   //-----------------------------
   unit_cell.atom[63].x   = 0.5;
   unit_cell.atom[63].y   = 0;
   unit_cell.atom[63].z   = 0;
   unit_cell.atom[63].mat = 1;
   unit_cell.atom[63].lc  = 63;
   unit_cell.atom[63].hc  = 0;
   unit_cell.atom[63].nm  = false;
   //-----------------------------
   unit_cell.atom[64].x   = 0.371;
   unit_cell.atom[64].y   = 0.629;
   unit_cell.atom[64].z   = 0;
   unit_cell.atom[64].mat = 2;
   unit_cell.atom[64].lc  = 64;
   unit_cell.atom[64].hc  = 0;
   unit_cell.atom[64].nm  = true;
   //-----------------------------
   unit_cell.atom[65].x   = 0.629;
   unit_cell.atom[65].y   = 0.371;
   unit_cell.atom[65].z   = 0;
   unit_cell.atom[65].mat = 2;
   unit_cell.atom[65].lc  = 65;
   unit_cell.atom[65].hc  = 0;
   unit_cell.atom[65].nm  = true;
   //-----------------------------
   unit_cell.atom[66].x   = 0.871;
   unit_cell.atom[66].y   = 0.871;
   unit_cell.atom[66].z   = 0.5;
   unit_cell.atom[66].mat = 2;
   unit_cell.atom[66].lc  = 66;
   unit_cell.atom[66].hc  = 2;
   unit_cell.atom[66].nm  = true;
   //-----------------------------
   unit_cell.atom[67].x   = 0.129;
   unit_cell.atom[67].y   = 0.129;
   unit_cell.atom[67].z   = 0.5;
   unit_cell.atom[67].mat = 2;
   unit_cell.atom[67].lc  = 67;
   unit_cell.atom[67].hc  = 2;
   unit_cell.atom[67].nm  = true;

   return;

}

// This function replicates a unit cell to a 3x3x3 grid

std::vector<local::atom_t> replicate(const unitcell::unit_cell_t& unit_cell){
   const double ucsx = unit_cell.dimensions[0];
   const double ucsy = unit_cell.dimensions[1];
   const double ucsz = unit_cell.dimensions[2];

   std::vector<local::atom_t> ratoms;

   for(int x = 0; x < 3; ++x){
      for(int y = 0; y < 3; ++y){
         for(int z = 0; z < 3; ++z){
            for(int a=0; a < unit_cell.atom.size(); ++a){
               local::atom_t tmp;
               tmp.mat = unit_cell.atom[a].mat;
               tmp.lc = unit_cell.atom[a].lc;
               tmp.hc = unit_cell.atom[a].hc;
               tmp.x = (unit_cell.atom[a].x + double(x))*ucsx;
               tmp.y = (unit_cell.atom[a].y + double(y))*ucsy;
               tmp.z = (unit_cell.atom[a].z + double(z))*ucsz;
               tmp.id = a;
               tmp.idx = x; // unit cell id
               tmp.idy = y;
               tmp.idz = z;
               tmp.nm = unit_cell.atom[a].nm;
               ratoms.push_back(tmp);
            }
         }
      }
   }
   return ratoms;
}

// Takes in 3x3x3 unit cell grid and outputs all interactions for central unit
// cell atoms
// Note all arguments should be in unit cell side length units
std::vector<local::interaction_t> radial_Neighbour_Bin_Assignment(
   std::vector<local::atom_t> ratoms, double radial_increment){
   
   // Define data structure to be output
   std::vector<local::interaction_t> interaction_data;
   
   // Find central unit cell magnetic atoms
   for(int i = 0; i < ratoms.size(); ++i){
      if ( (ratoms[i].idx == 1) && (ratoms[i].idy == 1) && (ratoms[i].idz == 1) &&
      !(ratoms[i].nm) ){
         // Ensuring other atom of interaction is magnetic and j!=i, and that only Fe-Fe and Fe-Nd interactions are considered.
         for(int j = 0; j < ratoms.size(); ++j){
            if (!(ratoms[j].nm) && j!=i && !(ratoms[i].mat == ratoms[j].mat && ratoms[i].mat == 0)){
               // temp interaction_t to assign and append data to output
               local::interaction_t tmp;
               // calculate interatomic radius and subsequently the radial bin
               const double rx = ratoms[j].x - ratoms[i].x;
               const double ry = ratoms[j].y - ratoms[i].y;
               const double rz = ratoms[j].z - ratoms[i].z;
               double inter_radius = sqrt(rx*rx + ry*ry + rz*rz);
               // Assign values to output vector.
               tmp.bin = (int) std::floor(inter_radius / radial_increment);
               tmp.mat = ratoms[i].mat;
               tmp.lc = ratoms[i].lc;
               tmp.hc = ratoms[i].hc;
               tmp.inter_mat = ratoms[j].mat;
               tmp.inter_hc = ratoms[j].hc;
               interaction_data.push_back(tmp);
            }
         }
      }
   }
   return interaction_data;
}

// These output data to a file so that a histogram can be plotted (not used)
void histogram_Data_Generator_General ( std::vector<local::interaction_t> interactions, double max_radius , double radial_bin) {
   int max_bin = std::floor( ( max_radius + radial_bin / 2 ) / radial_bin );
   std::ofstream file;
   file.open ( "General_Histogram.txt" );
   for ( int bin = 0; bin < max_bin; ++bin ) {
      int bin_count = 0;
      for ( int i = 0; i < interactions.size(); ++i ) { 
         if ( interactions[i].bin == bin ) {
            bin_count++;
         }
      }
      // Divide by no. of atoms to get per atom count.
      double per_Atom_Bin_Count = (double) bin_count / (double) 64;
      // Output bin no. and bin count. to file or create histogram object to plot from C++.
      file << (double) bin - 0.5 << " " << per_Atom_Bin_Count << std::endl;
   }
   file.close();
}

void histogram_Data_Generator_NdFe ( std::vector<local::interaction_t> interactions, double max_radius , double radial_bin) {
   int max_bin = std::floor( ( max_radius + radial_bin / 2 ) / radial_bin );
   std::ofstream file;
   file.open ( "NdFe_Histogram.txt" );
   for ( int bin = 0; bin < max_bin; ++bin ) {
      int bin_count = 0;
      for ( int i = 0; i < interactions.size(); ++i ) { 
         if ( interactions[i].bin == bin && interactions[i].mat == 0 && interactions[i].inter_mat == 1) {
            bin_count++;
         }
      }
      // Divide bin count by no. of Nd atoms to get per atom number.
      double per_Atom_Bin_Count = (double) bin_count / (double) 8;
      // Output bin no. and bin count. to file or create histogram object to plot from C++.
      file << (double) bin - 0.5 << " " << per_Atom_Bin_Count << std::endl;
   }
   file.close();
}

void histogram_Data_Generator_HC0 ( std::vector<local::interaction_t> interactions, double max_radius , double radial_bin) {
   int max_bin = std::floor( ( max_radius + radial_bin / 2 ) / radial_bin );
   std::ofstream file;
   file.open ( "HC0_Histogram.txt" );
   for ( int bin = 0; bin < max_bin; ++bin ) {
      int bin_count = 0;
      for ( int i = 0; i < interactions.size(); ++i ) { 
         if ( interactions[i].bin == bin && interactions[i].hc == 0 && (interactions[i].inter_hc == 1 || interactions[i].inter_hc == 3)) {
            bin_count++;
         }
      }
      // Divide bin count by no. of HC0 atoms to get per atom number.
      double per_Atom_Bin_Count = (double) bin_count / (double) 8;
      // Output bin no. and bin count. to file or create histogram object to plot from C++.
      file << (double) bin - 0.5 << " " << per_Atom_Bin_Count << std::endl;
   }
   file.close();
}

void histogram_Data_Generator_HC2 ( std::vector<local::interaction_t> interactions, double max_radius , double radial_bin) {
   int max_bin = std::floor( ( max_radius + radial_bin / 2 ) / radial_bin );
   std::ofstream file;
   file.open ( "HC2_Histogram.txt" );
   for ( int bin = 0; bin < max_bin; ++bin ) {
      int bin_count = 0;
      for ( int i = 0; i < interactions.size(); ++i ) { 
         if ( interactions[i].bin == bin && interactions[i].hc == 2 && (interactions[i].inter_hc == 1 || interactions[i].inter_hc == 3)) {
            bin_count++;
         }
      }
      // Divide bin count by no. of HC2 atoms to get per atom number.
      double per_Atom_Bin_Count = (double) bin_count / (double) 8;
      // Output bin no. and bin count. to file or create histogram object to plot from C++.
      file << (double) bin - 0.5 << " " << per_Atom_Bin_Count << std::endl;
   }
   file.close();
}

// rdf functions

void rdf_General (std::vector<local::interaction_t> interactions, double max_radius, double radial_bin) {
   int max_bin = std::floor( ( max_radius + radial_bin / 2 ) / radial_bin );
   std::ofstream file;
   file.open ( "rdf_General.txt" );
   //int tot_count = 0;
   for ( int bin = 0; bin < max_bin; ++bin ) {
      int bin_count = 0;
      for ( int i = 0; i < interactions.size(); ++i ) { 
         if ( interactions[i].bin == bin ) {
            bin_count++;
            //tot_count++;
         }
      }
      // Divide by no. of atoms to get per atom count.
      double rdf_Value = (double) bin_count / (double) 64;
      // Normalise for standard rdf.
      // Output radial distribution function data for gnuplot.
      file << std::fixed;
      file << std::setprecision(3);
      file << ((double) bin + 0.5)* 8.8 * radial_bin;
      file << " " << rdf_Value << std::endl;
   }
   file.close();
   //std::cout << (double) tot_count / 64.0 << std::endl;
}

void rdf_Nd (std::vector<local::interaction_t> interactions, double max_radius , double radial_bin) {
   int max_bin = std::floor( ( max_radius + radial_bin / 2 ) / radial_bin );
   std::ofstream file;
   file.open ( "rdf_Nd.txt" );
   //int tot_count = 0;
   for ( int bin = 0; bin < max_bin; ++bin ) {
      int bin_count = 0;
      for ( int i = 0; i < interactions.size(); ++i ) { 
         if ( interactions[i].bin == bin && interactions[i].mat == 0 && interactions[i].inter_mat == 1) {
            bin_count++;
            //tot_count++;
         }
      }
      // Divide bin count by no. of Nd atoms to get per atom value
      double rdf_Value = (double) bin_count / (double) 8;
      // Output radial distribution function data for gnuplot
      file << std::fixed;
      file << std::setprecision(3);
      file << ((double) bin + 0.5) * 8.8 * radial_bin;
      file << " " << rdf_Value << std::endl;
   }
   file.close();
   //std::cout << (double) tot_count / 8.0 << std::endl;
}

void rdf_HC0 ( std::vector<local::interaction_t> interactions, double max_radius , double radial_bin) {
   int max_bin = std::floor( ( max_radius + radial_bin / 2 ) / radial_bin );
   std::ofstream file;
   //int tot_count = 0;
   file.open ( "rdf_HC0.txt" );
   for ( int bin = 0; bin < max_bin; ++bin ) {
      int bin_count = 0;
      for ( int i = 0; i < interactions.size(); ++i ) { 
         if ( interactions[i].bin == bin && interactions[i].hc == 0 && (interactions[i].inter_hc == 1 || interactions[i].inter_hc == 3)) {
            bin_count++;
            //tot_count++;
         }
      }
      // Divide bin count by no. of HC0 atoms to get per atom number
      double rdf_Value = (double) bin_count / (double) 8;
      // Output radial distribution function data for gnuplot
      file << std::fixed;
      file << std::setprecision(3);
      file << ((double) bin + 0.5) * 8.8 * radial_bin;
      file << " " << rdf_Value << std::endl;
   }
   file.close();
   //std::cout << (double) tot_count / 8.0 << std::endl;
}

void rdf_lc0 (std::vector<local::interaction_t> interactions, double max_radius , double radial_bin) {
   int max_bin = std::floor( ( max_radius + radial_bin / 2 ) / radial_bin );
   std::ofstream file;
   file.open ( "rdf_lc0.txt" );
   //int tot_count = 0;
   for ( int bin = 0; bin < max_bin; ++bin ) {
      int bin_count = 0;
      for ( int i = 0; i < interactions.size(); ++i ) { 
         if ( interactions[i].bin == bin && interactions[i].lc == 0 && interactions[i].inter_mat == 1) {
            bin_count++;
            //tot_count++;
         }
      }
      // Divide bin count by no. of Nd atoms to get per atom value
      double rdf_Value = (double) bin_count;
      // Output radial distribution function data for gnuplot
      file << std::fixed;
      file << std::setprecision(3);
      file << ((double) bin + 0.5) * 8.8 * radial_bin;
      file << " " << rdf_Value << std::endl;
   }
   file.close();
   //std::cout << (double) tot_count / 8.0 << std::endl;
}

void rdf_lc60 (std::vector<local::interaction_t> interactions, double max_radius , double radial_bin) {
   int max_bin = std::floor( ( max_radius + radial_bin / 2 ) / radial_bin );
   std::ofstream file;
   file.open ( "rdf_lc60.txt" );
   //int tot_count = 0;
   for ( int bin = 0; bin < max_bin; ++bin ) {
      int bin_count = 0;
      for ( int i = 0; i < interactions.size(); ++i ) { 
         if ( interactions[i].bin == bin && interactions[i].lc == 60 && interactions[i].inter_mat == 1) {
            bin_count++;
            //tot_count++;
         }
      }
      // Divide bin count by no. of Nd atoms to get per atom value
      double rdf_Value = (double) bin_count;
      // Output radial distribution function data for gnuplot
      file << std::fixed;
      file << std::setprecision(3);
      file << ((double) bin + 0.5) * 8.8 * radial_bin;
      file << " " << rdf_Value << std::endl;
   }
   file.close();
   //std::cout << (double) tot_count / 8.0 << std::endl;
}

// Function to find the least number of nearest neighbours an atom in the bulk can have.

void bulkLeastNearestNeighbours( double radial_cutoff, double bin_size, std::vector<local::interaction_t> interactions ) {
   int lowest_Count = 10000;
   int radial_cutoff_bin = std::floor((radial_cutoff + bin_size / 2) / bin_size);
   for ( int i = 0; i < 64; ++i ){
      int count = 0;
      for ( int j = 0; j < interactions.size(); ++j ) {
         if ( interactions[j].lc == i && interactions[j].bin < radial_cutoff_bin) {
            count++;
         }
      }
      if ( count < lowest_Count ) {
         lowest_Count = count;
      }
   }
   std::cout << lowest_Count << std::endl;
}

//------------------------------------------------------------------------------
// Run the code
//------------------------------------------------------------------------------

int main(){
   // This function considers atoms in the central unit cell and creates a radial
   // probability distribution function for nearest neighbours.

   //--------------------------------------------------------------------------
   // Create 3x3x3 unit cell structure
   //--------------------------------------------------------------------------

   // Initialise a NdFeB unit cell structure
   unitcell::unit_cell_t NdFeB;
   build_NdFeB_Cell(NdFeB);

   // Replicate NdFeB to 3x3x3 grid
   std::vector<local::atom_t> ratoms = replicate(NdFeB);

   //-----------------------------------------------------------------------------
   // Create and plot radial probability distribution functions
   //-----------------------------------------------------------------------------

   // Set histogram radial binning width
   double radial_bin = 0.0001;

   // Calculate radial bin data
   std::vector<local::interaction_t> interactions = radial_Neighbour_Bin_Assignment(ratoms, radial_bin);
   
   // Output data

   //for ( int i = 0; i < interactions.size(); ++i ) {
   //   std::cout << interactions[i].lc << " " << interactions[i].mat << " " << interactions[i].hc << " " << interactions[i].inter_mat << " " << interactions[i].inter_hc << " " << interactions[i].bin << std::endl;
   //}

   // Output rdf data

//   histogram_Data_Generator_General(interactions, 0.5, radial_bin );

//   histogram_Data_Generator_NdFe(interactions, 0.5, radial_bin );

//   histogram_Data_Generator_HC0(interactions, 0.5, radial_bin );

//   histogram_Data_Generator_HC2(interactions, 0.5, radial_bin );

   bulkLeastNearestNeighbours(0.4725, radial_bin, interactions );

   //rdf_General(interactions, 1.0, radial_bin);

   //rdf_Nd(interactions, 1.0, radial_bin);

   //rdf_HC0(interactions, 1.0, radial_bin);

   //rdf_lc0(interactions, 1.0, radial_bin);

   //rdf_lc60(interactions, 1.0, radial_bin);

   return 0;
}