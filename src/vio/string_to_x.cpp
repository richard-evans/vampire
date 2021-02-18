//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2018. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <sstream>

// Vampire headers
#include "vio.hpp"

// vio module headers
#include "internal.hpp"

// vin namespace
namespace vin{

//------------------------------------------------------------
// Simple function to convert string to unit64_t
//------------------------------------------------------------
uint64_t str_to_uint64(std::string input_str){

   // intermediate double variable
   double value_dbl = 0.0;

   // load value into std::sstream for safe type conversion
   std::stringstream value_ss(input_str);

   // read value into double
   value_ss >> value_dbl;

   // truncate negative numbers to zero
   if(value_dbl < 0.0){
      value_dbl = 0.0;
   }

   // cast double to uint64_t
   uint64_t value = uint64_t(value_dbl);

   // return unit64_t
   return value;

}

//------------------------------------------------------------------------------
// Simple function to convert string to double
//------------------------------------------------------------------------------
double str_to_double(std::string input_str){

   // intermediate double variable
   double value = 0.0;

   // load value into std::sstream for safe type conversion
   std::stringstream value_ss(input_str);

   // read value into double
   value_ss >> value;

   // return double value
   return value;

}

} // end of namespace vin
