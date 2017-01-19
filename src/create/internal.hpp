#ifndef CREATE_INTERNAL_H_
#define CREATE_INTERNAL_H_
//-----------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2016. All rights reserved.
//
//-----------------------------------------------------------------------------

//---------------------------------------------------------------------
// Defines shared internal data structures and functions for the
// create implementation. These functions should
// not be accessed outside of the create module.
//---------------------------------------------------------------------

// transitional arrangement - need old create header for catom_t
// Should eventually move class definition to this file
#include "create.hpp"

// Vampire headers
#include "material.hpp"
#include "mtrand.hpp"

namespace create{
   namespace internal{

      //-----------------------------------------------------------------------------
      // Internal data type definitions
      //-----------------------------------------------------------------------------

      // alloy datatypes
      enum host_alloy_d_t { homogeneous, random, granular };
      enum slave_alloy_d_t { native, reciprocal, uniform };

      // simple class for slave material properties
      class slave_material_t{

      public:
         double fraction;
         double variance;
         slave_alloy_d_t slave_alloy_distribution;
         // constructor
         slave_material_t ():
         	fraction(0.0),
            variance(0.1),
            slave_alloy_distribution(native)
            {};
      };

      //-----------------------------------------------------------------------------
      // materials class for storing create material parameters
      //-----------------------------------------------------------------------------
      class mp_t{

      private:

      public:
         // variables
         bool alloy_master; // flag specifying if material is host material
         double host_alloy_smoothness;
         double host_alloy_scale;
         bool save_host_alloy_profile;
         std::string save_file_name;
         bool oxygen;
         bool oxygen_terminated;
         double oxygen_coordination_range;
         host_alloy_d_t host_alloy_distribution; // enum specifying type of alloy distribution
         std::vector<slave_material_t> slave_material; // array of slave alloys for host
         // constructor
         mp_t ():
         	alloy_master(false),
            host_alloy_smoothness(2.0),
            host_alloy_scale (50.0),
            save_host_alloy_profile(false),
            save_file_name(""),
            oxygen(false),
            oxygen_terminated(false),
            oxygen_coordination_range(1.0), // Unit cells
            host_alloy_distribution(internal::homogeneous)
            {
               // resize array of slave materials
               slave_material.resize(mp::max_materials);
            };
      };

      //-----------------------------------------------------------------------------
      // Internal shared variables used for creation
      //-----------------------------------------------------------------------------
      extern std::vector<create::internal::mp_t> mp; // array of material properties
      extern MTRand grnd; // general random number generator for create functions

      extern double faceted_particle_100_radius; // 100 facet radius
      extern double faceted_particle_110_radius; // 110 facet radius
      extern double faceted_particle_111_radius; // 111 facet radius
      extern double cone_angle; // angle of cone to truncate cylinder

      extern bool calculate_oxygen_termination; // determines whether oxygen termination should be calculated

      //-----------------------------------------------------------------------------
      // Internal functions for create module
      //-----------------------------------------------------------------------------
      extern void alloy(std::vector<cs::catom_t> & catom_array);
      extern void faceted(std::vector<double>& particle_origin, std::vector<cs::catom_t> & catom_array, const int grain);
      extern void clear_atoms(std::vector<cs::catom_t> &);
      extern void cone(std::vector<double>& particle_origin, std::vector<cs::catom_t> & catom_array, const int grain);
      extern void oxygen(std::vector<cs::catom_t> & catom_array);

   } // end of internal namespace
} // end of create namespace

#endif //CREATE_INTERNAL_H_
