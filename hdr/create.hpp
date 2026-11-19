//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2009-2012. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
#ifndef CREATE_H_
#define CREATE_H_
///
/// @file
/// @brief Contains the cs namespace header.
///
/// @details This is the detailed description of the funtion of this file
///
/// @section notes Implementation Notes
/// This is a list of other notes, not related to functionality but rather to implementation.
/// Also include references, formulae and other notes here.
///
///
/// @section info File Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    05/03/2010
/// @internal
///   Created:    05/03/2010
///   Revision:     ---
///=====================================================================================
///
#include <list>
#include <string>
#include <vector>
#include <cmath>

// Vampire headers
#include "create_atoms_class.hpp" // class definition for atoms in create module
#include "neighbours.hpp"
#include "unitcell.hpp"

/// @namespace
/// @brief Contains all functions and data associated with system creation in vampire.
///
/// @internal
///=====================================================================================
///

namespace cs{

   // System Dimensions
   extern double system_dimensions[3];
   extern bool pbc[3];
   extern bool SelectMaterialByGeometry;
   extern unsigned int total_num_unit_cells[3];
   extern unsigned int local_num_unit_cells[3];

   // System Parameters
   extern int particle_creation_parity;
   extern double particle_scale;
   extern double particle_spacing;
   extern double particle_array_offset_x; /// Offset particle array along x-direction;
   extern double particle_array_offset_y; /// Offset particle array along y-direction;
   extern double particle_shape_factor_x; /// Normalised particle shape
   extern double particle_shape_factor_y; /// Normalised particle shape
   extern double particle_shape_factor_z; /// Normalised particle shape

   // Other directives and flags
   extern bool single_spin;
   extern int system_creation_flags[10];
   extern bool fill_core_shell;
   extern bool core_shell_particles;

   // Variables for interfacial roughness control
   extern bool interfacial_roughness;
   extern bool interfacial_roughness_local_height_field;
   extern int interfacial_roughness_type; /// Sets peaks (1), troughs (-1) or both (0)
   extern unsigned int interfacial_roughness_random_seed;
   extern unsigned int interfacial_roughness_seed_count; /// Number of seeds
   extern double interfacial_roughness_height_field_resolution; // Angstroms
   extern double interfacial_roughness_mean_seed_radius; // Angstroms
   extern double interfacial_roughness_seed_radius_variance; // Variance as fraction of mean radius
   extern double interfacial_roughness_mean_seed_height; // Angstroms
   extern double interfacial_roughness_seed_height_max; // Angstroms

   // Variables for multilayer system
   extern bool multilayers;
   extern int num_multilayers;
   extern bool multilayer_height_category; // enable height categorization by multilayer number

   class neighbour_t {
   public:

      int nn; // atom id of neighbour
      int i; // interaction type of neighbour

      double vx; // vector between atoms i->j
      double vy;
      double vz;

   };

   extern uc::unit_cell_t unit_cell;

   // Structure for storing non-magnetic atom data
   struct nm_atom_t{
      double x;
      double y;
      double z;
      int mat;
      int cat;
      int grain;
   };

  // Array for storing non-magnetic atoms
  extern std::vector<nm_atom_t> non_magnetic_atoms_array;

/// @brief This is the brief (one line only) description of the function.
///
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    05/03/2010
///
/// @param[in] input variable
/// @param[out] ouput variable
/// @param[in,out] input/output variable
/// @return variable returned from the function
///
/// @internal
///   Created:    05/03/2010
///   Revision:     ---
///=====================================================================================
///
int create();

/// @brief This is the brief (one line only) description of the function.
///
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    05/03/2010
///
/// @param[in] input variable
/// @param[out] ouput variable
/// @param[in,out] input/output variable
/// @return variable returned from the function
///
/// @internal
///   Created:    05/03/2010
///   Revision:     ---
///=====================================================================================
///
int create_crystal_structure(std::vector<cs::catom_t> &);

int voronoi_film(std::vector<cs::catom_t> &);

void generate_multilayers(std::vector<cs::catom_t> & catom_array);

}

//------------------------------------------------------------------------------
// new create module functions
//------------------------------------------------------------------------------
namespace create{

   // Variable for total number of atoms that are not filler
   extern int num_total_atoms_non_filler;

   // Key of each atom (core, boundary and halo): built from its unit cell (supercell
   // coordinates, wrapped for periodic images) and its number within the unit cell, so it
   // does not depend on how the system is decomposed. Valid only if atom_keys_valid.
   extern std::vector<uint64_t> atom_keys;
   extern bool atom_keys_valid;
   void set_atom_keys(const std::vector<cs::catom_t>& catom_array);


   // Functions
   void initialize();
   int create_system_type(std::vector<cs::catom_t> &);
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);
   double get_material_height_min(const int material);
   double get_material_height_max(const int material);
   bool get_material_sub_fill(const int material);
   int get_material_unit_cell_category(const int material);
   double get_material_substructure_nucleation_height(const int material);

   // A material's core-shell radius (mp::material[mat].core_shell_size) and
   // the materials list ordered by that radius, largest first, so that
   // iterating in order and overwriting on each match paints the smallest,
   // innermost core last, on top of the larger shells enclosing it. Used by
   // any driver (particle shapes, granular films) that supports core-shell
   // materials.
   struct core_radius_t{ int mat; double radius; };
   bool compare_radius(core_radius_t first, core_radius_t second);
   std::list<core_radius_t> sorted_core_shell_materials();

   // Sorts catom_array in place by grain number, for improved cache
   // locality of later per-grain processing. Used by every system-creation
   // path that assigns atoms to grains (particle arrays, the nanoparticle
   // agglomerate, granular films).
   void sort_atoms_by_grain(std::vector<cs::catom_t> & catom_array);

} // end of namespace create

#endif /*CREATE_H_*/
