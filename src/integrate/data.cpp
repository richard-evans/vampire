//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <vector>

// Vampire headers
#include "integrate.hpp"

// integrate module headers
#include "internal.hpp"

namespace integrate{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   integrator_t integrator = integrator_t::llg_heun; // variable to specify integrator

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside integrate module
      //------------------------------------------------------------------------

      bool enabled; // bool to enable module

      std::vector<internal::mp_t> mp; // array of material properties

      std::vector<double> lsf_second_order_coefficient; // LSF Landau coefficients per material (J)
      std::vector<double> lsf_fourth_order_coefficient;
      std::vector<double> lsf_sixth_order_coefficient;

   } // end of internal namespace

   //------------------------------------------------------------------------------
   // Function to get LSF Landau coefficients (J) for a material
   //------------------------------------------------------------------------------
   void get_lsf_coefficients(const int material, double& A, double& B, double& C){
      A = internal::lsf_second_order_coefficient[material];
      B = internal::lsf_fourth_order_coefficient[material];
      C = internal::lsf_sixth_order_coefficient[material];
      return;
   }

} // end of integrate namespace

