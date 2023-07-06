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
      int init_num_materials = internal::mp.size();

      // if no anisotropy constants initialised, then make sure anisotropy array is the correct size
      if(init_num_materials == 0) internal::mp.resize(mu_s_array.size());

      // set actual number of materials
      const int num_materials = internal::mp.size();

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
      // Fourth order biaxial (simple version)
      if(internal::enable_biaxial_fourth_order_simple){
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
      // Fourth order rotational
      if(internal::enable_fourth_order_rotational){
         internal::k4r.resize(num_materials);
         for(int m = 0; m < num_materials; m++) internal::k4r[m] = internal::mp[m].k4r * inverse_mu_s[m];
      }

      //------------------------------------------------------------------------
      // 4th order triaxial
      //------------------------------------------------------------------------
      if(internal::enable_triaxial_fourth_order){

         internal::ku4_triaxial_vector_x.resize(num_materials);
         internal::ku4_triaxial_vector_y.resize(num_materials);
         internal::ku4_triaxial_vector_z.resize(num_materials);

         for(int m = 0; m < num_materials; m++) {
           internal::ku4_triaxial_vector_x[m] = internal::ku4_triaxial_vector_x[m] * inverse_mu_s[m];
           internal::ku4_triaxial_vector_y[m] = internal::ku4_triaxial_vector_y[m] * inverse_mu_s[m];
           internal::ku4_triaxial_vector_z[m] = internal::ku4_triaxial_vector_z[m] * inverse_mu_s[m];
         }

         //check orthogonality for fourth order basis sets
         for (int mat = 0; mat < num_materials; mat ++ ){
            if (!internal::triaxial_fourth_order_fixed_basis[mat]){
               double e1[3] = {internal::ku4_triaxial_basis1x[mat],internal::ku4_triaxial_basis1y[mat],internal::ku4_triaxial_basis1z[mat]};
               double e2[3] = {internal::ku4_triaxial_basis2x[mat],internal::ku4_triaxial_basis2y[mat],internal::ku4_triaxial_basis2z[mat]};
               double e3[3] = {internal::ku4_triaxial_basis3x[mat],internal::ku4_triaxial_basis3y[mat],internal::ku4_triaxial_basis3z[mat]};
               double cross[3] = {0,0,0};
               double onedottwo, onedotthree;

               //if basis 1 has been set in the material file
               if (e1[0] != 0 || e1[1] != 0 || e1[2] != 0){
                  //if basis 2 has been set in the material file
                  if (e2[0] != 0 || e2[1] != 0 || e2[2] != 0){
                     //are basis 1 and 2 orthogonal (1.2 = 0?)
                     onedottwo = e1[0]*e2[0] + e1[1]*e2[1] + e1[2]*e2[2];
                     if (onedottwo < 0.05){
                        //if so work out the 3rd orthoogonal basis set
                        cross[0] = e1[1] * e2[2] - e1[2] * e2[1];
                        cross[1] = e1[2] * e2[0] - e1[0] * e2[2];
                        cross[2] = e1[0] * e2[1] - e1[1] * e2[0];
                        //was the third already set in the input file?
                        if ((e3[0] != 0 || e3[1] != 0 || e3[2] != 0)){
                           //std::cout << "SET" << "\t" << e3[0] << '\t' << e3[1] << '\t' << e3[2] << std::endl;
                           //does it equal the one calcuated if not set to be the orthogonal basis set and print to logfile that i ahve done that
                           if ((cross[0] - e3[0] < 0.05 && cross[1] - e3[1] < 0.05 && cross[2] - e3[2] < 0.05)){}
                           else{
                              e3[0] = cross[0];
                              e3[1] = cross[1];
                              e3[2] = cross[2];
                              std::cerr << "Basis 3 is not orthogonal to basis 1,2 in material " << mat << " changing basis 3 to: (" <<  cross[0] << "," <<cross[1] << "," <<cross[2] << ")" << std::endl;
                              zlog << zTs() << "Basis 3 is not orthogonal to basis 1,2 in material " << mat << " changing basis 3 to: (" <<  cross[0] << "," <<cross[1] << "," <<cross[2] << ")" << std::endl;
                           }
                        }
                        else {
                           e3[0] = cross[0];
                           e3[1] = cross[1];
                           e3[2] = cross[2];
                           //std::cerr << "Basis 3 for material  " << mat << " is set to: (" <<  cross[0] << "," <<cross[1] << "," <<cross[2] << ")" << std::endl;
                           zlog << zTs() << "Basis 3 for material  " << mat << " is set to: (" <<  cross[0] << "," <<cross[1] << "," <<cross[2] << ")" << std::endl;
                        }
                     }
                     else {
                        std::cerr << "Basis 1,2 are not orthogonal for material:" << mat << std::endl;
                        zlog << zTs() << "Basis 1,2 are not orthogonal for material:" << mat << std::endl;
                        err::vexit();
                     }

                  }
                  else if ((e3[0] != 0 || e3[1] != 0 || e3[2] != 0)){
                     onedotthree = e1[0]*e3[0] + e1[1]*e3[1] + e1[2]*e3[2];
                     if (onedotthree < 0.05){
                        //if so work out the 3rd orthoogonal basis set
                        cross[0] = e1[1] * e3[2] - e1[2] * e3[1];
                        cross[1] = e1[2] * e3[0] - e1[0] * e3[2];
                        cross[2] = e1[0] * e3[1] - e1[1] * e3[0];
                        e2[0] = cross[0];
                        e2[1] = cross[1];
                        e2[2] = cross[2];
                        //std::cerr << "Basis 2 for material  " << mat << " is set to: (" <<  cross[0] << "," <<cross[1] << "," <<cross[2] << ")" << std::endl;
                        zlog << zTs() << "Basis 2 for material  " << mat << " is set to: (" <<  cross[0] << "," <<cross[1] << "," <<cross[2] << ")" << std::endl;

                     }
                     else {
                        std::cerr << "Basis 1,3 are not orthogonal for material:" << mat << std::endl;
                        zlog << zTs() << "Basis 1,3 are not orthogonal for material:" << mat << std::endl;
                        err::vexit();
                     }

                  }
                  else {
                     std::cerr << "Only one basis vector set for material:" << mat << " Please specify another basis vector" << std::endl;
                     zlog << zTs() << "Only one basis vector set for material:" << mat << " Please specify another basis vector" << std::endl;
                     err::vexit();
                  }
               }

            }
            else{
               internal::ku4_triaxial_basis1x[mat] = 1;
               internal::ku4_triaxial_basis1y[mat] = 0;
               internal::ku4_triaxial_basis1z[mat] = 0;
               internal::ku4_triaxial_basis2x[mat] = 0;
               internal::ku4_triaxial_basis2y[mat] = 1;
               internal::ku4_triaxial_basis2z[mat] = 0;
               internal::ku4_triaxial_basis3x[mat] = 0;
               internal::ku4_triaxial_basis3y[mat] = 0;
               internal::ku4_triaxial_basis3z[mat] = 1;

            }
         }

         //---------------------------------------------------------------------
         // override which version of triaxial anisotropy is needed
         //---------------------------------------------------------------------
         bool fixed_basis_fourth_order = true;
         for (int mat = 0; mat < num_materials; mat ++ ){
            if (internal::triaxial_fourth_order_fixed_basis[mat] == false){
               fixed_basis_fourth_order = false;
            }
         }
         // if any material requires rotated basis set, then use fancy function
         if(fixed_basis_fourth_order == false){
            internal::enable_triaxial_fourth_order = false;
            internal::enable_triaxial_fourth_order_rotated = true;
         }

      }

      //------------------------------------------------------------------------
      // Second order triaxial
      //------------------------------------------------------------------------
      if(internal::enable_triaxial_anisotropy){

         internal::ku_triaxial_vector_x.resize(num_materials);
         internal::ku_triaxial_vector_y.resize(num_materials);
         internal::ku_triaxial_vector_z.resize(num_materials);

         for(int m = 0; m < num_materials; m++) {
            internal::ku_triaxial_vector_x[m] = internal::ku_triaxial_vector_x[m] * inverse_mu_s[m];
            internal::ku_triaxial_vector_y[m] = internal::ku_triaxial_vector_y[m] * inverse_mu_s[m];
            internal::ku_triaxial_vector_z[m] = internal::ku_triaxial_vector_z[m] * inverse_mu_s[m];
         }

         //check orthogonality for second order basis sets
         for (int mat = 0; mat < num_materials; mat ++ ){
            if (!internal::triaxial_second_order_fixed_basis[mat]){
               double e1[3] = {internal::ku_triaxial_basis1x[mat],internal::ku_triaxial_basis1y[mat],internal::ku_triaxial_basis1z[mat]};
               double e2[3] = {internal::ku_triaxial_basis2x[mat],internal::ku_triaxial_basis2y[mat],internal::ku_triaxial_basis2z[mat]};
               double e3[3] = {internal::ku_triaxial_basis3x[mat],internal::ku_triaxial_basis3y[mat],internal::ku_triaxial_basis3z[mat]};
               double cross[3] = {0,0,0};
               double onedottwo, onedotthree;

               //if basis 1 has been set in the material file
               if (e1[0] != 0 || e1[1] != 0 || e1[2] != 0){
                  //if basis 2 has been set in the material file
                  if (e2[0] != 0 || e2[1] != 0 || e2[2] != 0){
                     //are basis 1 and 2 orthogonal (1.2 = 0?)
                     onedottwo = e1[0]*e2[0] + e1[1]*e2[1] + e1[2]*e2[2];
                     if (onedottwo < 0.05){
                        //if so work out the 3rd orthoogonal basis set
                        cross[0] = e1[1] * e2[2] - e1[2] * e2[1];
                        cross[1] = e1[2] * e2[0] - e1[0] * e2[2];
                        cross[2] = e1[0] * e2[1] - e1[1] * e2[0];
                        //was the third already set in the input file?
                        if ((e3[0] != 0 || e3[1] != 0 || e3[2] != 0)){
                           //std::cout << "SET" << "\t" << e3[0] << '\t' << e3[1] << '\t' << e3[2] << std::endl;
                           //does it equal the one calcuated if not set to be the orthogonal basis set and print to logfile that i ahve done that
                           if ((cross[0] - e3[0] < 0.05 && cross[1] - e3[1] < 0.05 && cross[2] - e3[2] < 0.05)){}
                           else{
                              e3[0] = cross[0];
                              e3[1] = cross[1];
                              e3[2] = cross[2];
                              std::cerr << "Basis 3 is not orthogonal to basis 1,2 in material " << mat << " changing basis 3 to: (" <<  cross[0] << "," <<cross[1] << "," <<cross[2] << ")" << std::endl;
                              zlog << zTs() << "Basis 3 is not orthogonal to basis 1,2 in material " << mat << " changing basis 3 to: (" <<  cross[0] << "," <<cross[1] << "," <<cross[2] << ")" << std::endl;
                           }
                        }
                        else {
                           e3[0] = cross[0];
                           e3[1] = cross[1];
                           e3[2] = cross[2];
                           //std::cerr << "Basis 3 for material  " << mat << " is set to: (" <<  cross[0] << "," <<cross[1] << "," <<cross[2] << ")" << std::endl;
                           zlog << zTs() << "Basis 3 for material  " << mat << " is set to: (" <<  cross[0] << "," <<cross[1] << "," <<cross[2] << ")" << std::endl;
                        }
                     }
                     else {
                        std::cerr << "Basis 1,2 are not orthogonal for material:" << mat << std::endl;
                        zlog << zTs() << "Basis 1,2 are not orthogonal for material:" << mat << std::endl;
                        err::vexit();
                     }

                  }
                  else if ((e3[0] != 0 || e3[1] != 0 || e3[2] != 0)){
                     onedotthree = e1[0]*e3[0] + e1[1]*e3[1] + e1[2]*e3[2];
                     if (onedotthree < 0.05){
                        //if so work out the 3rd orthoogonal basis set
                        cross[0] = e1[1] * e3[2] - e1[2] * e3[1];
                        cross[1] = e1[2] * e3[0] - e1[0] * e3[2];
                        cross[2] = e1[0] * e3[1] - e1[1] * e3[0];
                        e2[0] = cross[0];
                        e2[1] = cross[1];
                        e2[2] = cross[2];
                        std::cerr << "Basis 2 for material  " << mat << " is set to: (" <<  cross[0] << "," <<cross[1] << "," <<cross[2] << ")" << std::endl;
                        zlog << zTs() << "Basis 2 for material  " << mat << " is set to: (" <<  cross[0] << "," <<cross[1] << "," <<cross[2] << ")" << std::endl;

                     }
                     else {
                        std::cerr << "Basis 1,3 are not orthogonal for material:" << mat << std::endl;
                        zlog << zTs() << "Basis 1,3 are not orthogonal for material:" << mat << std::endl;
                        err::vexit();
                     }

                  }
                  else {
                     std::cerr << "Only one basis vector set for material:" << mat << " Please specify another basis vector" << std::endl;
                     zlog << zTs() << "Only one basis vector set for material:" << mat << " Please specify another basis vector" << std::endl;
                     err::vexit();
                  }
               }

            }
            else{
               internal::ku_triaxial_basis1x[mat] = 1;
               internal::ku_triaxial_basis1y[mat] = 0;
               internal::ku_triaxial_basis1z[mat] = 0;
               internal::ku_triaxial_basis2x[mat] = 0;
               internal::ku_triaxial_basis2y[mat] = 1;
               internal::ku_triaxial_basis2z[mat] = 0;
               internal::ku_triaxial_basis3x[mat] = 0;
               internal::ku_triaxial_basis3y[mat] = 0;
               internal::ku_triaxial_basis3z[mat] = 1;

            }
         }

         //---------------------------------------------------------------------
         // override which version of triaxial anisotropy is needed
         //---------------------------------------------------------------------
         bool fixed_basis_second_order = true;
         for (int mat = 0; mat < num_materials; mat ++ ){
            if (!internal::triaxial_second_order_fixed_basis[mat]){
               fixed_basis_second_order = false;
            }
         }
         // if any material requires rotated basis set, then use fancy function
         if(fixed_basis_second_order == false){
            internal::enable_triaxial_anisotropy = false;
            internal::enable_triaxial_anisotropy_rotated = true;
         }

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
         internal::klattice.resize(num_materials);
         internal::klattice_array.resize(num_materials); // anisoptropy constant

         // loop over all materials and set up lattice anisotropy constants
         for(int m = 0; m < num_materials; m++){

            // save anisotropy constant to unrolled array in units of tesla
            internal::klattice[m] = internal::mp[m].k_lattice * inverse_mu_s[m];

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
