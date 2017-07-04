//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef EXCHANGE_H_
#define EXCHANGE_H_

// C++ standard library headers
#include <string>

// Vampire headers
#include "create.hpp"
#include "exchange.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions for exchange module
//--------------------------------------------------------------------------------
namespace exchange{

   //-----------------------------------------------------------------------------
   // Function to initialise exchange module
   //-----------------------------------------------------------------------------
   void initialize(std::vector<std::vector <cs::neighbour_t> >& cneighbourlist);

   //---------------------------------------------------------------------------
   // Function to process input file parameters for exchange module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

} // end of exchange namespace

#endif //EXCHANGE_H_
