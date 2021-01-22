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


      //---------------------------------------------------------------------
      // get number of materials for simulation
      //---------------------------------------------------------------------
      unsigned int init_num_materials = internal::mp.size();

      // if no anisotropy constants initialised, then make sure anisotropy array is the correct size
      if(init_num_materials == 0) internal::mp.resize(mu_s_array.size());

      // set actual number of materials
      const unsigned int num_materials = internal::mp.size();

      // output informative message
      zlog << zTs() << "Initialising data structures for anisotropy calculation for " << num_materials << " materials" << std::endl;

      // check for prior initialisation
      if (internal::initialised){
         zlog << zTs() << "Warning: Anisotropy calculation already initialised. Continuing." << std::endl;
         return;
      }

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
      if(internal::enable_cubic_fourth_order || internal::enable_cubic_fourth_order_rotation){
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

      for(int m = 0; m < num_materials; m++){

         // unroll uniaxial easy axes
         internal::ku_vector[m].x = internal::mp[m].ku_vector[0];
         internal::ku_vector[m].y = internal::mp[m].ku_vector[1];
         internal::ku_vector[m].z = internal::mp[m].ku_vector[2];

      }

      //---------------------------------------------------------------------
      // initialise rotated axis directions for each material
      //---------------------------------------------------------------------

      for(int mat = 0; mat < num_materials; mat++){

         // Vectors defining the easy axis in cubic anisotropy (Roberto was here)
         double e1[3] = { internal::mp[mat].kc_vector1[0],
                          internal::mp[mat].kc_vector1[1],
                          internal::mp[mat].kc_vector1[2] };

         double e2[3] = { internal::mp[mat].kc_vector2[0],
                          internal::mp[mat].kc_vector2[1],
                          internal::mp[mat].kc_vector2[2] };

         // calculate e3 as vector product e1 ^ e2
         double e3[3] = { (internal::mp[mat].kc_vector1[1]*internal::mp[mat].kc_vector2[2] - internal::mp[mat].kc_vector1[2]*internal::mp[mat].kc_vector2[1]),
                          (internal::mp[mat].kc_vector1[2]*internal::mp[mat].kc_vector2[0] - internal::mp[mat].kc_vector1[0]*internal::mp[mat].kc_vector2[2]),
                          (internal::mp[mat].kc_vector1[0]*internal::mp[mat].kc_vector2[1] - internal::mp[mat].kc_vector1[1]*internal::mp[mat].kc_vector2[0])};

         // Calculate vector lengths
         double mod_e1 = sqrt(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]);
         double mod_e2 = sqrt(e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2]);
         double mod_e3 = sqrt(e3[0]*e3[0] + e3[1]*e3[1] + e3[2]*e3[2]);

         // check for zero vectors and exit with error
         if(mod_e1 < 1e-9 || mod_e2 < 1e-9 || mod_e3 < 1e-9){
            std::cerr << "Error! Rotated cubic anisotropy vectors for material " << mat << " are not orthogonal. Exiting" << std::endl;
            zlog << zTs() << "Error! Rotated cubic anisotropy vectors for material " << mat << " are not orthogonal. Exiting" << std::endl;
            err::vexit();
         }

         // normalise vectors to unit length
         internal::mp[mat].kc_vector1[0] = e1[0] / mod_e1;
         internal::mp[mat].kc_vector1[1] = e1[1] / mod_e1;
         internal::mp[mat].kc_vector1[2] = e1[2] / mod_e1;

         internal::mp[mat].kc_vector2[0] = e2[0] / mod_e2;
         internal::mp[mat].kc_vector2[1] = e2[1] / mod_e2;
         internal::mp[mat].kc_vector2[2] = e2[2] / mod_e2;

         internal::mp[mat].kc_vector3[0] = e3[0] / mod_e3;
         internal::mp[mat].kc_vector3[1] = e3[1] / mod_e3;
         internal::mp[mat].kc_vector3[2] = e3[2] / mod_e3;

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
