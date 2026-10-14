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

   // constraint angles for constrained Monte Carlo (sim:constraint-angle-* keywords)
   bool   constraint_rotation=false; /// enables rotation of spins to new constraint direction
   bool   constraint_phi_changed=false; /// flag to note change in phi
   double constraint_phi=0.0; /// Constrained minimisation vector (azimuthal) [degrees]
   double constraint_phi_min=0.0; /// loop angle min [degrees]
   double constraint_phi_max=0.0; /// loop angle max [degrees]
   double constraint_phi_delta=5.0; /// loop angle delta [degrees]

   bool   constraint_theta_changed=false;
   double constraint_theta=0.0; /// Constrained minimisation vector (rotational) [degrees]
   double constraint_theta_min=0.0; /// loop angle min [degrees]
   double constraint_theta_max=0.0; // loop angle max [degrees]
   double constraint_theta_delta=5.0; /// loop angle delta [degrees]

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
      std::vector<double> Sold(3);
      std::vector<double> Snew(3);

      //MC-MPI variables
      std::vector<std::vector<int> > c_octants; //Core atoms of each octant
      std::vector<std::vector<int> > b_octants; //Boundary atoms of each octant


   } // end of internal namespace

} // end of montecarlo namespace
