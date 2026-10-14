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
#include "integrate.hpp"
#include "material.hpp"
#include "montecarlo.hpp"
#include "sim.hpp"
#include "errors.hpp"
#include "vio.hpp"

// montecarlo module headers
#include "internal.hpp"

namespace montecarlo{

   //----------------------------------------------------------------------------
   // Function to initialize montecarlo module
   //----------------------------------------------------------------------------
   void initialize(const int num_atoms, const int num_grains, const std::vector<int>& grain_array){
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

      //------------------------------------------------------------------------
      // Set up masked cmc with grain constraints
      //------------------------------------------------------------------------
      if(cmc::constrain_by_grain){

         // first check for hybrid cmc
         if(integrate::integrator == integrate::integrator_t::hybrid_cmc){
            // implementation not yet done - give message to user
            err::fatal("Error - grain constraints with materials are not yet supported with cmc. Please use sim:integrator = contrained-monte-carlo. Exiting");
         }

         // standard mode with global constraint (integrator = cmc)
         if(integrate::integrator == integrate::integrator_t::cmc){
            // these should be refactored into the mc module
            const double constraint_phi   = montecarlo::constraint_phi;
            const double constraint_theta = montecarlo::constraint_theta;

            // deep copy the grain numbers to mask
            std::vector<int> constraint_mask = grain_array;

            // all atoms are constrained
            std::vector<bool> constrained(num_atoms, true);

            // generate the same phi, theta for all atoms
            std::vector<double> phi_theta_constraints(2*num_grains);

            for(int i = 0; i < num_grains; i++){
               phi_theta_constraints[2*i+0] = constraint_phi;
               phi_theta_constraints[2*i+1] = constraint_theta;
            }

            // initialise masked cmc
            montecarlo::initialise_masked_cmc_mc(num_grains, constraint_mask, constrained, phi_theta_constraints);

            // now reset integrator to hybrid
            integrate::integrator = integrate::integrator_t::hybrid_cmc;

            // and turn on masked cmc integrator
            cmc::masked_cmc = true;

         }

      } // end of grain constraint check

      return;

   }

} // end of montecarlo namespace
