//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sam Westmoreland and Richard Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "anisotropy.hpp"

// anisotropy module headers
#include "internal.hpp"

namespace anisotropy{

   //--------------------------------------------------------------------------------
   // Function to get second order uniaxial anisotropy constant for a given material
   //--------------------------------------------------------------------------------
   double get_anisotropy_constant(const int material){
      return internal::mp[material].ku2;
   }

   //--------------------------------------------------------------------------------
   // Function to get second order uniaxial anisotropy constant for a given material
   //--------------------------------------------------------------------------------
   double get_ku2(const int material){
      return internal::mp[material].ku2;
   }

   //--------------------------------------------------------------------------------
   // Function to get fourth order uniaxial anisotropy constant for a given material
   //--------------------------------------------------------------------------------
   double get_ku4(const int material){
      return internal::mp[material].ku4;
   }

   //--------------------------------------------------------------------------------
   // Function to get sixth order uniaxial anisotropy constant for a given material
   //--------------------------------------------------------------------------------
   double get_ku6(const int material){
      return internal::mp[material].ku6;
   }

   //--------------------------------------------------------------------------------
   // Function to get fourth order cubic anisotropy constant for a given material
   //--------------------------------------------------------------------------------
   double get_kc4(const int material){
      return internal::mp[material].kc4;
   }

   //--------------------------------------------------------------------------------
   // Function to get sixth order cubic anisotropy constant for a given material
   //--------------------------------------------------------------------------------
   double get_kc6(const int material){
      return internal::mp[material].kc6;
   }

   //--------------------------------------------------------------------------------
   // Function to get unit vector defining axis for uniaxial anisotropy for a material
   //--------------------------------------------------------------------------------
   std::vector<double> get_ku_vector(const int material){
      return internal::mp[material].ku_vector;
   }


} // end of anisotropy namespace
