//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sam Westmoreland and Richard Evans 2017. All rights reserved.
//
//   Email: sw766@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "anisotropy.hpp"

// anisotropy module headers
#include "internal.hpp"

namespace anisotropy{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside anisotropy module
      //------------------------------------------------------------------------

      std::vector<internal::mp_t> mp; // array of material properties

      bool initialised = false; // variable to determine if module has been initialised

      bool enable_second_order_tensor = false; // Flag to enable calculation of second order tensor anisotropy
      bool enable_fourth_order_tensor = false; // Flag to enable calculation of second order tensor anisotropy
      bool enable_sixth_order_tensor  = false; // Flag to enable calculation of second order tensor anisotropy

      bool enable_neel_anisotropy     = false; // Flag to turn on Neel anisotropy calculation (memory intensive at startup)
      bool enable_lattice_anisotropy  = false; // Flag to turn on lattice anisotropy calculation
      bool enable_random_anisotropy   = false; // Flag to enable random anisitropy initialisation

      // arrays for storing 1D collapsed tensors
      std::vector<double> second_order_tensor(0);
      std::vector<double> fourth_order_tensor(0);
      std::vector<double> sixth_order_tensor(0);

   	//bool identify_surface_atoms = false; // flag to identify surface atoms in config coordinate file
      bool native_neel_anisotropy_threshold  = false; // enables site-dependent surface threshold
   	unsigned int neel_anisotropy_threshold = 123456789; // global threshold for surface atoms
      double nearest_neighbour_distance      = 1.e9; // Control surface anisotropy nearest neighbour distance

      // arrays for storing unrolled parameters for lattice anisotropy
      std::vector<double> klattice_array(0); // anisoptropy constant
      std::vector<evec_t> elattice_array(0); // easy axis

   } // end of internal namespace

} // end of anisotropy namespace
