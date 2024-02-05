//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Jack Collings, Sam Westmoreland and Richard Evans 2017.
//   All rights reserved.
//
//   Email: sw766@york.ac.uk
//          jbc525@york.ac.uk
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

      std::vector<internal::mp_t> mp(0); // array of material properties

      bool initialised = false; // variable to determine if module has been initialised

      bool enable_neel_anisotropy     = false; // Flag to turn on Neel anisotropy calculation (memory intensive at startup)
      bool enable_lattice_anisotropy  = false; // Flag to turn on lattice anisotropy calculation
      bool enable_random_anisotropy   = false; // Flag to enable random anisitropy initialisation

      // Flags for spherical-harmonic decomposition of magnetocrystalline anisotropy
      bool enable_uniaxial_second_order      = false; // Flag to enable calculation of second order anisotropy
      bool enable_rotational_2_1_order       = false; // Flag to enable calculation of second order theta first ordrer phi anisotropy
      bool enable_rotational_2_1_order_odd   = false; // Flag to enable calculation of second order theta first order phi odd anisotropy
      bool enable_rotational_2_2_order       = false; // Flag to enable calculation of second order theta second order phi anisotropy
      bool enable_rotational_2_2_order_odd   = false; // Flag to enable calculation of second order theta second order phi odd anisotropy
      bool enable_uniaxial_fourth_order      = false; // Flag to enable calculation of fourth order anisotropy
      bool enable_rotational_4_1_order_odd   = false; // Flag to enable calculation of fourth order theta first order phi odd anisotropy
      bool enable_rotational_4_2_order       = false; // Flag to enable calculation of fourth order theta second order phi anisotorpy
      bool enable_rotational_4_2_order_odd   = false; // Flag to enable calculation of fourth order theta second order phi odd anisotropy
      bool enable_rotational_4_3_order       = false; // Flag to enable calculation of fourth order theta third order phi anisotropy
      bool enable_rotational_4_3_order_odd   = false; // Flag to enable calculation of fourth order theta third order phi odd anisotropy
      bool enable_rotational_4_4_order       = false; // Flag to enable calculation of fourth order theta fourth order phi anisotropy
      bool enable_rotational_4_4_order_odd   = false; // Flag to enable calculation of fourth order theta fourth order phi odd anisotropy
      bool enable_uniaxial_sixth_order       = false; // Flag to enable calculation of sixth order anisotropy
      bool enable_rotational_6_2_order       = false; // Flag to enable calculation of sixth order theta second order phi anisotropy
      bool enable_rotational_6_2_order_odd   = false; // Flag to enable calculation of sixth order theta second order phi odd anisotropy
      bool enable_rotational_6_4_order       = false; // Flag to enable calculation of sixth order theta fourth order phi anisotropy
      bool enable_rotational_6_4_order_odd   = false; // Flag to enable calculation of sixth order theta fourth order phi odd anisotropy
      bool enable_rotational_6_6_order       = false; // Flag to enable calculation of sixth order theta sixth order phi anisotropy
      bool enable_rotational_6_6_order_odd   = false; // Flag to enable calculation of sixth order theta sixth order phi odd anisotropy

      bool enable_biaxial_fourth_order_simple = false; // Flag to enable calculation of the simple version of the fourth order anisotropy

      bool enable_cubic_fourth_order    = false; // Flag to enable calculation of fourth order cubic anisotropy
      bool enable_cubic_sixth_order     = false; // Flag to enable calculation of sixth order cubic  anisotropy
      bool enable_cubic_fourth_order_rotation = false; // Flag to enable calculation of rotated cubic anisotropy

      bool  enable_triaxial_anisotropy = false;
      bool  enable_triaxial_fourth_order = false;
      bool  enable_triaxial_anisotropy_rotated = false;
      bool  enable_triaxial_fourth_order_rotated = false;

      // array for storing 1D second order collapsed tensor for Neel anisotropy
      std::vector<double> neel_tensor(0);

      // arrays for storing unrolled anisotropy constants in Tesla
      std::vector<double> ku2(0);
      std::vector<double> k2r1(0);
      std::vector<double> k2r1_odd(0);
      std::vector<double> k2r2(0);
      std::vector<double> k2r2_odd(0);
      std::vector<double> ku4(0);
      std::vector<double> k4r1_odd(0);
      std::vector<double> k4r2(0);
      std::vector<double> k4r2_odd(0);
      std::vector<double> k4r3(0);
      std::vector<double> k4r3_odd(0);
      std::vector<double> k4r4(0);
      std::vector<double> k4r4_odd(0);
      std::vector<double> ku6(0);
      std::vector<double> k6r2(0);
      std::vector<double> k6r2_odd(0);
      std::vector<double> k6r4(0);
      std::vector<double> k6r4_odd(0);
      std::vector<double> k6r6(0);
      std::vector<double> k6r6_odd(0);

      std::vector<double> kc4(0);
      std::vector<double> kc6(0);

      // unrolled arrays for storing easy axes for each material
      std::vector<evec_t> ku_vector(0); // 001 easy axis direction
      std::vector<evec_t> kr_vector(0); // 100 rotational axis direction
      std::vector<evec_t> kl_vector(0); // 010 last axis direction

      std::vector<evec_t> u1_vector(0); // Unit vector along [110]
      std::vector<evec_t> u2_vector(0); // Unit vector along [1-10]

      std::vector<double> ku_triaxial_vector_x(100,0); // unit vector defining axis for triaxial anisotropy
      std::vector<double> ku_triaxial_vector_y(100,0); //
      std::vector<double> ku_triaxial_vector_z(100,0); //

      std::vector<double> ku4_triaxial_vector_x(100,0); // unit vector defining axis for triaxial anisotropy
      std::vector<double> ku4_triaxial_vector_y(100,0); //
      std::vector<double> ku4_triaxial_vector_z(100,0); //

      //basis vectors for second order triaxial - must be orthogonality
      std::vector < double > ku_triaxial_basis1x(100,0.0);
      std::vector < double > ku_triaxial_basis1y(100,0.0);
      std::vector < double > ku_triaxial_basis1z(100,0.0);

      std::vector < double > ku_triaxial_basis2x(100,0.0);
      std::vector < double > ku_triaxial_basis2y(100,0.0);
      std::vector < double > ku_triaxial_basis2z(100,0.0);

      std::vector < double > ku_triaxial_basis3x(100,0.0);
      std::vector < double > ku_triaxial_basis3y(100,0.0);
      std::vector < double > ku_triaxial_basis3z(100,0.0);

      //basis vectors for fourth order triaxial - must be orthogonality
      std::vector < double > ku4_triaxial_basis1x(100,0.0);
      std::vector < double > ku4_triaxial_basis1y(100,0.0);
      std::vector < double > ku4_triaxial_basis1z(100,0.0);

      std::vector < double > ku4_triaxial_basis2x(100,0.0);
      std::vector < double > ku4_triaxial_basis2y(100,0.0);
      std::vector < double > ku4_triaxial_basis2z(100,0.0);

      std::vector < double > ku4_triaxial_basis3x(100,0.0);
      std::vector < double > ku4_triaxial_basis3y(100,0.0);
      std::vector < double > ku4_triaxial_basis3z(100,0.0);

      std::vector<bool> triaxial_second_order_fixed_basis(100,true);
      std::vector<bool> triaxial_fourth_order_fixed_basis(100,true);

      bool native_neel_anisotropy_threshold  = false; // enables site-dependent surface threshold
   	unsigned int neel_anisotropy_threshold = 123456789; // global threshold for surface atoms
      double nearest_neighbour_distance      = 1.e9; // Control surface anisotropy nearest neighbour distance
      bool neel_range_dependent              = false; // Enable range dependent Neel anisotropy Lij = L0 exp(-F(r-r0)/r0 )
      double neel_exponential_range          = 2.5; // r0 value for range dependence of Neel anisotropy
      double neel_exponential_factor         = 5.53; // F value for range dependence of Neel anisotropy (default assumes nnn fraction of 10% of nn value)

      // arrays for storing unrolled parameters for lattice anisotropy
      std::vector<double> klattice_array(0); // anisoptropy constant

   } // end of internal namespace

} // end of anisotropy namespace
