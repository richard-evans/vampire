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

      // Initialise tensors for each atom
      internal::second_order_tensor.resize( 9 * num_atoms, 0.0 );
      internal::fourth_order_tensor.resize( 9 * num_atoms, 0.0 );
      internal::sixth_order_tensor.resize ( 9 * num_atoms, 0.0 );

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

         internal::neel_anisotropy(num_atoms);

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

      //---------------------------------------------------------------------
      // Populate sixth order tensor
      //---------------------------------------------------------------------
      if (internal::enable_sixth_order_tensor){

         // Add uniaxial sixth order anisotropy (Ku3)
         //internal::uniaxial_sixth_order(num_atoms, atom_material_array, inverse_mu_s);

         // Add cubic sixth order anisotropy (Kc2)
         //internal::cubic_sixth_order(num_atoms, atom_material_array, inverse_mu_s);

      }

      //---------------------------------------------------------------------
      // initialise lattice anisotropy for each material
      //---------------------------------------------------------------------
      if(internal::enable_lattice_anisotropy){

         // get number of materials for simulation
         const unsigned int num_materials = mu_s_array.size();

         // arrays for storing unrolled parameters for lattice anisotropy
         internal::klattice_array.resize(num_materials); // anisoptropy constant
         internal::elattice_array.resize(num_materials); // easy axis

         // loop over all materials and set up lattice anisotropy constants
         for(int m = 0; m < num_materials; m++){

            // set up interpolation between temperature points
            internal::mp[m].lattice_anisotropy.set_interpolation_table();

            // unroll easy axes for lattice anisotropy calculation
            internal::elattice_array[m].x = internal::mp[m].ku_vector[0];
            internal::elattice_array[m].y = internal::mp[m].ku_vector[1];
            internal::elattice_array[m].z = internal::mp[m].ku_vector[2];

            // output interpolated data to file
            //internal::mp[m].lattice_anisotropy.output_interpolated_function(mat);

         }

      }

      internal::initialised = true;

      return;
   }

} // end of anisotropy namespace
