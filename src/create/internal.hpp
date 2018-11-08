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

      struct core_radius_t{
         int mat;
         double radius;
      };

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
         host_alloy_d_t host_alloy_distribution; // enum specifying type of alloy distribution
         std::vector<slave_material_t> slave_material; // array of slave alloys for host
         bool sub_fill; // flag to determine if material fills voided space in substructure
         double voronoi_grain_substructure_nucleation_height; // value determines start point of nucleated grains
         int unit_cell_category; // association of material to unit cell id
         double min; // minimum material height
         double max; // maximum material height
         double max_width;
         double min_width;

         // constructor
         mp_t ():
         	alloy_master(false),
            host_alloy_smoothness(2.0),
            host_alloy_scale (50.0),
            save_host_alloy_profile(false),
            save_file_name(""),
            host_alloy_distribution(internal::homogeneous),
            sub_fill(false),
            voronoi_grain_substructure_nucleation_height(0.0),
            unit_cell_category(0),
            min(0.0),
            max(1.0),
            min_width(0.0),
            max_width(1.0)
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

      extern int alloy_seed;  // random seed to control alloying of atoms
      extern int grain_seed;  // random seed to control grain structure generation
      extern int dilute_seed; // random seed to control dilution of atoms
      extern int mixing_seed; // random seed to control intermixing of atoms

      extern double faceted_particle_100_radius; // 100 facet radius
      extern double faceted_particle_110_radius; // 110 facet radius
      extern double faceted_particle_111_radius; // 111 facet radius
      extern double cone_angle; // angle of cone to truncate cylinder

      extern double voronoi_grain_size;
      extern double voronoi_grain_spacing;

      extern double bubble_radius;
      extern double bubble_nucleation_height;

      extern bool generate_voronoi_substructure;
      extern double voronoi_grain_substructure_crystallization_radius;
      extern double voronoi_grain_substructure_overlap_factor;

      extern bool select_material_by_z_height;

      //-----------------------------------------------------------------------------
      // Internal functions for create module
      //-----------------------------------------------------------------------------
      void set_atom_vars(std::vector<cs::catom_t> &, neighbours::list_t& bilinear, neighbours::list_t& biquadratic);


      extern void alloy(std::vector<cs::catom_t> & catom_array);
      extern void layers(std::vector<cs::catom_t> & catom_array);
      extern void roughness(std::vector<cs::catom_t> & catom_array);
      extern void bubble(std::vector<double>& particle_origin, std::vector<cs::catom_t> & catom_array, const int grain);
      extern void bulk(std::vector<cs::catom_t> & catom_array);
      extern void cone(std::vector<double>& particle_origin, std::vector<cs::catom_t> & catom_array, const int grain);
      extern void cube(std::vector<double>& particle_origin, std::vector<cs::catom_t> & catom_array, const int grain);
      extern void cylinder(std::vector<double>& particle_origin, std::vector<cs::catom_t> & catom_array, const int grain);
      extern void ellipsoid(std::vector<double>& particle_origin, std::vector<cs::catom_t> & catom_array, const int grain);
      extern void faceted(std::vector<double>& particle_origin, std::vector<cs::catom_t> & catom_array, const int grain);
      extern void sphere(std::vector<double>& particle_origin, std::vector<cs::catom_t> & catom_array, const int grain);
      extern void teardrop(std::vector<double>& particle_origin, std::vector<cs::catom_t> & catom_array, const int grain);
      extern void truncated_octahedron(std::vector<double>& particle_origin, std::vector<cs::catom_t> & catom_array, const int grain);

      extern void voronoi_substructure(std::vector<cs::catom_t> & catom_array);

      void voronoi_grain_rounding(std::vector <std::vector <double> > & grain_coord_array,
                                  std::vector <std::vector <std::vector <double> > > &  grain_vertices_array);

      void populate_vertex_points(std::vector <std::vector <double> > & grain_coord_array,
                                  std::vector <std::vector <std::vector <double> > > &  grain_vertices_array,
                                  bool include_boundary_grains);

      extern bool compare_radius(core_radius_t first,core_radius_t second);

      // MPI functions
      extern void copy_halo_atoms(std::vector<cs::catom_t> & catom_array);
      extern void identify_mpi_boundary_atoms(std::vector<cs::catom_t> & catom_array, neighbours::list_t & cneighbourlist);
      extern void mark_non_interacting_halo(std::vector<cs::catom_t>& catom_array);
      extern void sort_atoms_by_mpi_type(std::vector<cs::catom_t> & catom_array, neighbours::list_t& bilinear, neighbours::list_t& biquadratic);
      extern void init_mpi_comms(std::vector<cs::catom_t> & catom_array);

   } // end of internal namespace
} // end of create namespace

#endif //CREATE_INTERNAL_H_
