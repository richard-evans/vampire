//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sam Westmoreland and Richard Evans 2017. All rights reserved.
//
//   Email: sw766@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <vector>

// Vampire headers
#include "anisotropy.hpp"
//#include "atoms.hpp"
#include "errors.hpp"
#include "vio.hpp"
//#include "material.hpp"

// anisotropy module headers
#include "internal.hpp"

namespace anisotropy{

   //----------------------------------------------------------------------------
   // function to initialize anisotropy module
   //----------------------------------------------------------------------------
   void initialize (const unsigned int   num_atoms, // number of atoms
                    std::vector<int>&    atom_material_array, // atoms::atom_type_array
                    std::vector<double>& mu_s_array // array of magnetic moments
                   ){

      /* output informative message */
      zlog << zTs() << "Initialising data structures for anisotropy calculation." << std::endl;

      /* check for prior initialisation */
      if (internal::initialised){
         zlog << zTs() << "Warning: Anisotropy calculation already initialised. Continuing." << std::endl;
         return;
      }

      // Initialise tensor variables for each atom
      if(internal::enable_second_order_tensor) internal::second_order_tensor.resize( 9 * num_atoms, 0.0 );
      if(internal::enable_fourth_order_tensor) internal::fourth_order_tensor.resize( 9 * num_atoms, 0.0 );
      if(internal::enable_sixth_order_tensor)  internal::sixth_order_tensor.resize ( 9 * num_atoms, 0.0 );

      // Unroll inverse mu_S array for materials to convert Joules to Tesla
      const double mu_B = 9.27400915e-24; // Bohr magneton
      std::vector <double> inverse_mu_s(mu_s_array.size()); // array storing inverse spin moment in J/T
      for(int m = 0; m < mu_s_array.size(); m++) inverse_mu_s[m] = 1.0 / ( mu_s_array[m] * mu_B );

      //---------------------------------------------------------------------
      // Populate second order tensor
      //---------------------------------------------------------------------
      if (internal::enable_second_order_tensor){

         // Add uniaxial second order anisotropy (Ku1)
         internal::uniaxial_second_order(num_atoms, atom_material_array, inverse_mu_s);

      }

      //---------------------------------------------------------------------
      // Populate fourth order tensor
      //---------------------------------------------------------------------
      if (internal::enable_fourth_order_tensor){

         // Add uniaxial fourth order anisotropy (Ku2)
         internal::uniaxial_fourth_order(num_atoms, atom_material_array, inverse_mu_s);

         // Add cubic fourth order anisotropy (Kc1)
         internal::cubic_fourth_order(num_atoms, atom_material_array, inverse_mu_s);

      }

      internal::initialised = true;

      return;
   }

} // end of anisotropy namespace
