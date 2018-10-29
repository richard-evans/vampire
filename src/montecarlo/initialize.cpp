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
#include "material.hpp"
#include "montecarlo.hpp"

// montecarlo module headers
#include "internal.hpp"

namespace montecarlo{

   //----------------------------------------------------------------------------
   // Function to initialize montecarlo module
   //----------------------------------------------------------------------------
   void initialize(){
      //Copy materials data into internal namespace
      internal::num_materials = mp::num_materials;
      internal::temperature_rescaling_alpha.resize(mp::num_materials);
      internal::temperature_rescaling_Tc.resize(mp::num_materials);
      internal::mu_s_SI.resize(mp::num_materials);
      for(int i=0; i < mp::num_materials; i++) {
         internal::temperature_rescaling_alpha[i] = mp::material[i].temperature_rescaling_alpha;
         internal::temperature_rescaling_Tc[i] = mp::material[i].temperature_rescaling_Tc;
         internal::mu_s_SI[i] = mp::material[i].mu_s_SI;
      }

      //Initialize parallel mc variables
      mc_parallel_initialized = false;
      internal::c_octants.resize(8);
      internal::b_octants.resize(8);


      return;

   }

} // end of montecarlo namespace
