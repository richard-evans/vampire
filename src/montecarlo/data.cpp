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
#include <vector>

// Vampire headers
#include "montecarlo.hpp"

// montecarlo module headers
#include "internal.hpp"

namespace montecarlo{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   bool mc_parallel_initialized = false;
   bool lsf_mc_parallel_initialized = false;

   namespace cmc{
      bool masked_cmc = false;       // determine if generic masked constraint is required
      bool constrain_by_grain = false; // constrains by grain rather than globally
   }

   // Monte Carlo update algorithm
   algorithm_t algorithm = adaptive;

   namespace internal{
      //------------------------------------------------------------------------
      // Shared variables inside montecarlo module
      //------------------------------------------------------------------------

      // Materials variables
      int num_materials;
      std::vector<double> temperature_rescaling_alpha;
      std::vector<double> temperature_rescaling_Tc;
      std::vector<double> mu_s_SI;

      // MC Variables
      double delta_angle = 0.1;     // Tuned angle for Monte Carlo trial move
      double adaptive_sigma = 60.0; // sigma trial width for adaptive move
      //std::vector<double> adaptive_sigma_l = {0.1, 0.1};
      std::vector<double> Sold(3);
      std::vector<double> Snew(3);

      //MC-MPI variables
      std::vector<std::vector<int> > c_octants; //Core atoms of each octant
      std::vector<std::vector<int> > b_octants; //Boundary atoms of each octant


   } // end of internal namespace

} // end of montecarlo namespace
