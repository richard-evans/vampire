//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Adam Laverack and Richard Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "montecarlo.hpp"

// montecarlo module headers
#include "internal.hpp"

namespace montecarlo{

   //----------------------------------------------------------------------------
   // Function to initialize montecarlo module
   //----------------------------------------------------------------------------
   void initialize(){

      internal::mc_parallel_initialized = false;
      internal::c_octants.resize(8);
      internal::b_octants.resize(8);


      return;

   }

} // end of montecarlo namespace
