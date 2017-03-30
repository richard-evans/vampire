//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2016. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "create.hpp"

// Internal create header
#include "internal.hpp"

namespace create{

   //----------------------------------------------------------------------------
   // Shared variables used with main vampire code
   //---------------------------------------------------------------------------

      namespace internal{

         //----------------------------------------------------------------------------
         // Shared variables used within create module
         //---------------------------------------------------------------------------
         std::vector<create::internal::mp_t> mp; // array of material properties
         MTRand grnd; // general random number generator for create functions

         double faceted_particle_100_radius = 1.0; // 100 facet particle radius
         double faceted_particle_110_radius = 1.0; // 110 facet particle radius
         double faceted_particle_111_radius = 1.0; // 111 facet particle radius
         double cone_angle = 10.0; // factor to truncate cone

         double voronoi_grain_size = 50.0;
         double voronoi_grain_spacing = 10.0;

         double bubble_radius = 0.3333;
         double bubble_nucleation_height = 0.0;

         bool generate_voronoi_substructure = false;
         double voronoi_grain_substructure_crystallization_radius = 1.2;

      } // end of internal namespace

} // end of create namespace
