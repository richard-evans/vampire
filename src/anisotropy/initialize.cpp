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
#include "errors.hpp"
#include "vio.hpp"

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

      //---------------------------------------------------------------------
      // get number of materials for simulation
      //---------------------------------------------------------------------
      const unsigned int num_materials = mu_s_array.size();

      //---------------------------------------------------------------------
      // Unroll inverse mu_S array for materials to convert Joules to Tesla
      //---------------------------------------------------------------------
      const double mu_B = 9.27400915e-24; // Bohr magneton
      std::vector <double> inverse_mu_s(num_materials); // array storing inverse spin moment in J/T
      for(int m = 0; m < num_materials; m++) inverse_mu_s[m] = 1.0 / ( mu_s_array[m] * mu_B );

      //---------------------------------------------------------------------
      // Unroll material constants into arrays
      //---------------------------------------------------------------------
      // Second order uniaxial
      if(internal::enable_uniaxial_second_order){
         internal::ku2.resize(num_materials);
         for(int m = 0; m < num_materials; m++) internal::ku2[m] = internal::mp[m].ku2 * inverse_mu_s[m];
      }
      // Fourth order uniaxial
      if(internal::enable_uniaxial_fourth_order){
         internal::ku4.resize(num_materials);
         for(int m = 0; m < num_materials; m++) internal::ku4[m] = internal::mp[m].ku4 * inverse_mu_s[m];
      }
      // Sixth order uniaxial
      if(internal::enable_uniaxial_sixth_order){
         internal::ku6.resize(num_materials);
         for(int m = 0; m < num_materials; m++) internal::ku6[m] = internal::mp[m].ku6 * inverse_mu_s[m];
      }
      // Fourth order cubic
      if(internal::enable_cubic_fourth_order){
         internal::kc4.resize(num_materials);
         for(int m = 0; m < num_materials; m++) internal::kc4[m] = internal::mp[m].kc4 * inverse_mu_s[m];
      }
      // Sixth order cubic
      if(internal::enable_cubic_sixth_order){
         internal::kc6.resize(num_materials);
         for(int m = 0; m < num_materials; m++) internal::kc6[m] = internal::mp[m].kc6 * inverse_mu_s[m];
      }

      //---------------------------------------------------------------------
      // initialise axes for each material
      //---------------------------------------------------------------------
      internal::ku_vector.resize(num_materials);
      internal::kc_vector.resize(num_materials);

      for(int m = 0; m < num_materials; m++){

         // unroll uniaxial easy axes
         internal::ku_vector[m].x = internal::mp[m].ku_vector[0];
         internal::ku_vector[m].y = internal::mp[m].ku_vector[1];
         internal::ku_vector[m].z = internal::mp[m].ku_vector[2];

         // unroll uniaxial easy axes
         internal::kc_vector[m].x = internal::mp[m].kc_vector[0];
         internal::kc_vector[m].y = internal::mp[m].kc_vector[1];
         internal::kc_vector[m].z = internal::mp[m].kc_vector[2];

      }

      //---------------------------------------------------------------------
      // initialise lattice anisotropy for each material
      //---------------------------------------------------------------------
      if(internal::enable_lattice_anisotropy){

         // arrays for storing unrolled parameters for lattice anisotropy
         internal::klattice_array.resize(num_materials); // anisoptropy constant

         // loop over all materials and set up lattice anisotropy constants
         for(int m = 0; m < num_materials; m++){

            // set up interpolation between temperature points
            internal::mp[m].lattice_anisotropy.set_interpolation_table();

            // output interpolated data to file
            //internal::mp[m].lattice_anisotropy.output_interpolated_function(mat);

         }

      }

      //---------------------------------------------------------------------
      // set flag after initialization
      //---------------------------------------------------------------------
      internal::initialised = true;

      return;
   }

} // end of anisotropy namespace
