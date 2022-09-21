//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2022. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <fstream>
#include <iomanip>
#include <iostream>


// Vampire headers
#include "unitcell.hpp"

// unitcell module headers
#include "internal.hpp"

/*class atom_t {
public:
   double x; /// atom x-coordinate
   double y; /// atom y-coordinate
   double z; /// atom z-coordinate
   unsigned int mat; /// material
   unsigned int lc; /// lattice category
   unsigned int hc; /// height category
   unsigned int ni; /// number of interactions
   bool nm; // non-magnetic atom (no interactions are calculated)

   // constructor
   atom_t():
      x(0.0),
      y(0.0),
      z(0.0),
      mat(0),
      lc(0),
      hc(0),
      ni(0),
      nm(false)
   {
   };

};

//---------------------------------------------------------------------------
// Unit cell interaction class definition
//---------------------------------------------------------------------------
class interaction_t {
public:
   unsigned int i; /// atom unit cell id
   unsigned int j; /// neighbour atom unit cell id
   unsigned int mat_i; /// atom material category
   unsigned int mat_j; /// neighbour material category
   unsigned int shell; // shell number of interaction
   int dx; /// delta x in unit cells
   int dy; /// delta y in unit cells
   int dz; /// delta z in unit cells
   double rij; // interaction range (unit cells)
   double Jij[3][3]; /// Exchange tensor
};
class unit_cell_t {
public:

   double dimensions[3];
   double shape[3][3];
   double cutoff_radius; // nearest neighbours
   unsigned int interaction_range; /// maximum range in unit cells

   unsigned int lcsize; /// number of local categories
   unsigned int hcsize; /// number of height categories
   unsigned int surface_threshold; /// threshold for surface atoms

   // list of atoms in each unit cell
   std::vector <unitcell::atom_t> atom;

   unitcell::exchange_template_t bilinear;
   unitcell::exchange_template_t biquadratic;
   //exchange_template_t fourspin_interaction; // tbc

};

exchange::exchange_t exchange_type; // exchange type to use in simulation
bool use_material_exchange_constants; // flag to enable material exchange parameters
int num_unit_cell_atoms; // number of atoms in unit cell

// list of interactions in each unit cell
std::vector <unitcell::interaction_t> interaction;

// list of number of interactions from template for each atom in unit cell
std::vector <int> ni;*/

//------------------------------------------------------------------------------
// Templated functions to make output fixed width
//------------------------------------------------------------------------------
template<typename T> std::string fw(T var){
   std::stringstream ss, ssfw; // standard and fixed width forms
   ss << var;
   ssfw << std::setw(10) << var;
   // check that variable fits in size
   if(ss.str().size() <= 10) return ssfw.str();
   else return ss.str();
}

template<typename T> std::string fw5(T var){
   std::stringstream ss, ssfw; // standard and fixed width forms
   ss << var;
   ssfw << std::setw(5) << var;
   // check that variable fits in size
   if(ss.str().size() <= 5) return ssfw.str();
   else return ss.str();
}

namespace unitcell{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   namespace internal{



      //------------------------------------------------------------------------------
      // Function to write unit cell file to disk
      //------------------------------------------------------------------------------

      void write_unit_cell_file(unit_cell_t & uc){

         std::ofstream ofile("unitcell.ucf");

         ofile << "#-----------------------------------------------------------------------------------------" << std::endl;
         ofile << "# Generated unit cell file for vampire for structure " << uc::internal::crystal_structure << std::endl;
         ofile << "#-----------------------------------------------------------------------------------------" << std::endl;
         ofile << "# Unit cell dimensions (Angstroms)" << std::endl;
         ofile << "#-----------------------------------------------------------------------------------------" << std::endl;
         ofile << fw(uc.dimensions[0]) << "\t" << fw(uc.dimensions[1]) << "\t" << fw(uc.dimensions[2]) << std::endl;
         ofile << "#-----------------------------------------------------------------------------------------" << std::endl;
         ofile << "# Unit cell vectors: " << std::endl;
         ofile << "#-----------------------------------------------------------------------------------------" << std::endl;
         ofile << uc.shape[0][0] << "\t" << uc.shape[0][1] << "\t" << uc.shape[0][2] << std::endl;
         ofile << uc.shape[1][0] << "\t" << uc.shape[1][1] << "\t" << uc.shape[1][2] << std::endl;
         ofile << uc.shape[2][0] << "\t" << uc.shape[2][1] << "\t" << uc.shape[2][2] << std::endl;
         ofile << "#-----------------------------------------------------------------------------------------" << std::endl;
         ofile << "# Total number of atoms; atom id, cx cy cz, material, lattice category, height category" << std::endl;
         ofile << "#-----------------------------------------------------------------------------------------" << std::endl;
         ofile << uc.atom.size() << std::endl;
         for(int a = 0 ; a < uc.atom.size() ; a++ ){
            ofile << fw(a) << "\t" << fw(uc.atom[a].x) << "\t" << fw(uc.atom[a].y) << "\t" << fw(uc.atom[a].z) << "\t" <<
                     fw5(uc.atom[a].mat) << "\t" << fw5(uc.atom[a].lc) << "\t" << fw5(uc.atom[a].hc) << std::endl;
         }
         ofile << "#-----------------------------------------------------------------------------------------" << std::endl;
         ofile << "# Number of exchange interactions, exchange type; id, i j dx dy dz Jxx Jxy.. " << std::endl;
         ofile << "#-----------------------------------------------------------------------------------------" << std::endl;

         // set exchange names
         std::vector<std::string> exchange_names = {"isotropic", "vectorial", "tensorial" };
         std::vector<std::string> norm_exchange_names = {"normalised-isotropic", "normalised-vectorial", "normalised-tensorial" };

         // ternery function to grab exchange type
         std::string extype = uc.bilinear.use_material_exchange_constants ? norm_exchange_names[uc.bilinear.exchange_type] : exchange_names[uc.bilinear.exchange_type];

         ofile << fw5(uc.bilinear.interaction.size()) << "\t" << fw(extype) << std::endl;
         // loop over all exchange interations
         int count = 0;
         for(auto uci : uc.bilinear.interaction){
            ofile << fw5(count) << "\t" << fw5(uci.i) << "\t" << fw5(uci.j) << "\t"
                  << fw5(uci.dx) << "\t" << fw5(uci.dy) << "\t" << fw5(uci.dz) << "\t";
            switch(uc.bilinear.exchange_type){
               case(exchange::isotropic) : ofile << fw(uci.Jij[0][0]) << std::endl; break;
               case(exchange::vectorial) : ofile << fw(uci.Jij[0][0]) << "\t" << fw(uci.Jij[1][1]) << "\t" << fw(uci.Jij[2][2]) << std::endl; break;
               case(exchange::tensorial) : ofile << fw(uci.Jij[0][0]) << "\t" << fw(uci.Jij[0][1]) << "\t" << fw(uci.Jij[0][2]) << "\t"
                                                 << fw(uci.Jij[1][0]) << "\t" << fw(uci.Jij[1][1]) << "\t" << fw(uci.Jij[1][2]) << "\t"
                                                 << fw(uci.Jij[2][0]) << "\t" << fw(uci.Jij[2][1]) << "\t" << fw(uci.Jij[2][2]) << std::endl; break;
            }
            count++;
         }

         ofile.close();

         return;

      }

   } // end of internal namespace

} // end of unitcell namespace
