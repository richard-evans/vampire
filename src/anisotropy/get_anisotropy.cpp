//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Jack B. Collings (2022), Sam Westmoreland and Richard Evans 2017.
//   All rights reserved.
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
   // Function to get second order theta first order phi anisotropy constant for a
   // given material
   //--------------------------------------------------------------------------------
   double get_k2r1(const int material){
      return internal::mp[material].k2r1;
   }

   //--------------------------------------------------------------------------------
   // Function to get second order theta first order phi odd anisotropy constant for
   // a given material
   //--------------------------------------------------------------------------------
   double get_k2r1_odd(const int material){
      return internal::mp[material].k2r1_odd;
   }

   //--------------------------------------------------------------------------------
   // Function to get second order theta second order phi anisotropy constant for a
   // given material
   //--------------------------------------------------------------------------------
   double get_k2r2(const int material){
      return internal::mp[material].k2r2;
   }

   //--------------------------------------------------------------------------------
   // Function to get second order theta second order phi odd anisotropy constant for
   // a given material
   //--------------------------------------------------------------------------------
   double get_k2r2_odd(const int material){
      return internal::mp[material].k2r2_odd;
   }

   //--------------------------------------------------------------------------------
   // Function to get fourth order uniaxial anisotropy constant for a given material
   //--------------------------------------------------------------------------------
   double get_ku4(const int material){
      return internal::mp[material].ku4;
   }

   //--------------------------------------------------------------------------------
   // Function to get fourth order theta first order phi odd anisotropy constant for
   // a given material
   //--------------------------------------------------------------------------------
   double get_k4r1_odd(const int material){
      return internal::mp[material].k4r1_odd;
   }

   //--------------------------------------------------------------------------------
   // Function to get fourth order theta second order phi anisotropy constant for a
   // given material
   //--------------------------------------------------------------------------------
   double get_k4r2(const int material){
      return internal::mp[material].k4r2;
   }

   //--------------------------------------------------------------------------------
   // Function to get fourth order theta second order phi odd anisotropy constant for
   // a given material
   //--------------------------------------------------------------------------------
   double get_k4r2_odd(const int material){
      return internal::mp[material].k4r2_odd;
   }

   //--------------------------------------------------------------------------------
   // Function to get fourth order theta fourth order phi anisotropy constant for a
   // given material
   //--------------------------------------------------------------------------------
   double get_k4r4(const int material){
      return internal::mp[material].k4r4;
   }

   //--------------------------------------------------------------------------------
   // Function to get fourth order theta third order phi anisotropy constant for a
   // given material
   //--------------------------------------------------------------------------------
   double get_k4r3(const int material){
      return internal::mp[material].k4r3;
   }

   //--------------------------------------------------------------------------------
   // Function to get fourth order theta third order phi odd anisotropy constant for
   // a given material
   //--------------------------------------------------------------------------------
   double get_k4r3_odd(const int material){
      return internal::mp[material].k4r3_odd;
   }

   //--------------------------------------------------------------------------------
   // Function to get fourth order theta fourth order phi odd anisotropy constant for
   // a given material
   //--------------------------------------------------------------------------------
   double get_k4r4_odd(const int material){
      return internal::mp[material].k4r4_odd;
   }

   //--------------------------------------------------------------------------------
   // Function to get sixth order uniaxial anisotropy constant for a given material
   //--------------------------------------------------------------------------------
   double get_ku6(const int material){
      return internal::mp[material].ku6;
   }

   //--------------------------------------------------------------------------------
   // Function to get sixth order theta second order phi anisotropy constant for a
   // given material
   //--------------------------------------------------------------------------------
   double get_k6r2(const int material){
      return internal::mp[material].k6r2;
   }

   //--------------------------------------------------------------------------------
   // Function to get sixth order theta second order phi odd anisotropy constant for
   // a given material
   //--------------------------------------------------------------------------------
   double get_k6r2_odd(const int material){
      return internal::mp[material].k6r2_odd;
   }

   //--------------------------------------------------------------------------------
   // Function to get sixth order theta fourth order phi anisotropy constant for a
   // given material
   //--------------------------------------------------------------------------------
   double get_k6r4(const int material){
      return internal::mp[material].k6r4;
   }

   //--------------------------------------------------------------------------------
   // Function to get sixth order theta fourth order phi odd anisotropy constant for
   // a given material
   //--------------------------------------------------------------------------------
   double get_k6r4_odd(const int material){
      return internal::mp[material].k6r4_odd;
   }

   //--------------------------------------------------------------------------------
   // Function to get sixth order theta sixth order phi anisotropy constant for a
   // given material
   //--------------------------------------------------------------------------------
   double get_k6r6(const int material){
      return internal::mp[material].k6r4;
   }

   //--------------------------------------------------------------------------------
   // Function to get sixth order theta sixth order phi odd anisotropy constant for a
   // given material
   //--------------------------------------------------------------------------------
   double get_k6r6_odd(const int material){
      return internal::mp[material].k6r4_odd;
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
   // Function to get unit vector defining axis for uniaxial anisotropy for a
   // material
   //--------------------------------------------------------------------------------
   std::vector<double> get_ku_vector(const int material){
      return internal::mp[material].ku_vector;
   }

   std::vector<double> get_kr_vector(const int material){
      return internal::mp[material].kr_vector;
   }

   std::vector<double> get_kl_vector(const int material){
      return internal::mp[material].kl_vector;
   }


} // end of anisotropy namespace
