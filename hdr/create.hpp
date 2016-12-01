//-----------------------------------------------------------------------------
//
//  Vampire - A code for atomistic simulation of magnetic materials
//
//  Copyright (C) 2009-2012 R.F.L.Evans
//
//  Email:richard.evans@york.ac.uk
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
//  General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software Foundation,
//  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
//
// ----------------------------------------------------------------------------
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
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    05/03/2010
/// @internal
///	Created:		05/03/2010
///	Revision:	  ---
///=====================================================================================
///
#include <string>
#include <vector>
#include <cmath>

/// @namespace
/// @brief Contains all functions and data associated with system creation in vampire.
///
/// @internal
///=====================================================================================
///

namespace cs{

	// System Dimensions
	extern double system_dimensions[3];
	extern double unit_cell_size[3];
	extern bool pbc[3];
	extern bool SelectMaterialByZHeight;
	extern bool SelectMaterialByGeometry;
	extern unsigned int total_num_unit_cells[3];
	extern unsigned int local_num_unit_cells[3];
	extern std::string crystal_structure;

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
	extern std::string unit_cell_file;
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

	class unit_cell_atom_t {
	public:
		double x; /// atom x-coordinate
		double y; /// atom y-coordinate
		double z; /// atom z-coordinate
		unsigned int mat; /// material
		unsigned int lc; /// lattice category
		unsigned int hc; /// height category
		unsigned int ni; /// number of interactions
	};

	class unit_cell_interaction_t {
	public:
		unsigned int i; /// atom unit cell id
		unsigned int j; /// neighbour atom unit cell id
		int dx; /// delta x in unit cells
		int dy; /// delta y in unit cells
		int dz; /// delta z in unit cells
		double Jij[3][3]; /// Exchange tensor
	};

	class unit_cell_t {
	public:

		double dimensions[3];
		double shape[3][3];

		unsigned int lcsize; /// number of local categories
		unsigned int hcsize; /// number of height categories
		unsigned int interaction_range; /// maximum range in unit cells
		unsigned int surface_threshold; /// threshold for surface atoms
		int exchange_type; /// -1=isotropic(local material), 0=isotropic, 1=vector, or 2=tensor

		// list of atoms in each unit cell
		std::vector <unit_cell_atom_t> atom;

		// list of interactions in each unit cell
		std::vector <unit_cell_interaction_t> interaction;

	};

	class neighbour_t {
	public:

		int nn; // atom id of neighbour
		int i; // interaction type of neighbour

      double vx; // vector between atoms i->j
      double vy;
      double vz;

	};

	extern cs::unit_cell_t unit_cell;

   // Structure for storing non-magnetic atom data
   struct nm_atom_t{
      double x;
      double y;
      double z;
      int mat;
      std::string element;
   };

  // Array for storing non-magnetic atoms
  extern std::vector<nm_atom_t> non_magnetic_atoms_array;

	class catom_t {
		public:

			// Coordinates
			double x;
			double y;
			double z;

			// Flags
			bool include;

			// Integers
			int material;
			unsigned int uc_id;
			int uc_category;
			int lh_category;
			int grain;
			int supercell;
			int mpi_type;
			int mpi_cpuid;
			int mpi_atom_number;
			int mpi_old_atom_number;
			int scx;
			int scy;
			int scz;

			catom_t():
				x(0.0),
				y(0.0),
				z(0.0),
				include(false),
				material(0),
				uc_id(0),
				uc_category(0),
				lh_category(0),
				grain(0),
				supercell(0),
				mpi_type(0),
				mpi_cpuid(0),
				mpi_atom_number(0),
				mpi_old_atom_number(0),
				scx(0),
				scy(0),
				scz(0)
			{};
};
/// @brief This is the brief (one line only) description of the function.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
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
///	Created:		05/03/2010
///	Revision:	  ---
///=====================================================================================
///
int create();

/// @brief This is the brief (one line only) description of the function.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
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
///	Created:		05/03/2010
///	Revision:	  ---
///=====================================================================================
///
int create_crystal_structure(std::vector<cs::catom_t> &);

/// @brief This is the brief (one line only) description of the function.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
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
///	Created:		05/03/2010
///	Revision:	  ---
///=====================================================================================
///
int create_system_type(std::vector<cs::catom_t> &);

/// @brief This is the brief (one line only) description of the function.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
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
///	Created:		05/03/2010
///	Revision:	  ---
///=====================================================================================
///
int create_neighbourlist(std::vector<cs::catom_t> &, std::vector<std::vector <neighbour_t> > &);

/// @brief This is the brief (one line only) description of the function.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
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
///	Created:		05/03/2010
///	Revision:	  ---
///=====================================================================================
///
int set_atom_vars(std::vector<cs::catom_t> &, std::vector<std::vector <neighbour_t> > &);

int voronoi_film(std::vector<cs::catom_t> &);

/// @brief This is the brief (one line only) description of the function.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
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
///	Created:		05/03/2010
///	Revision:	  ---
///=====================================================================================
///
int bulk(std::vector<cs::catom_t> &);

/// @brief This is the brief (one line only) description of the function.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
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
///	Created:		05/03/2010
///	Revision:	  ---
///=====================================================================================
///
int cube(double[], std::vector<cs::catom_t> &,const int);

/// @brief This is the brief (one line only) description of the function.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
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
///	Created:		05/03/2010
///	Revision:	  ---
///=====================================================================================
///
int sphere(double[], std::vector<cs::catom_t> &,const int);

extern void ellipsoid(double[], std::vector<cs::catom_t> &,const int);

/// @brief This is the brief (one line only) description of the function.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
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
///	Created:		05/03/2010
///	Revision:	  ---
///=====================================================================================
///
int cylinder(double[], std::vector<cs::catom_t> &,const int);

/// @brief This is the brief (one line only) description of the function.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
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
///	Created:		05/03/2010
///	Revision:	  ---
///=====================================================================================
///
int truncated_octahedron(double[], std::vector<cs::catom_t> &,const int);
int tear_drop(double[], std::vector<cs::catom_t> &,const int);

int sort_atoms_by_grain(std::vector<cs::catom_t> &);

void roughness(std::vector<cs::catom_t> &);
void generate_multilayers(std::vector<cs::catom_t> & catom_array);

  // unit cell initialisation function
  void unit_cell_set(cs::unit_cell_t &);

}

//------------------------------------------------------------------------------
// new create module functions
//------------------------------------------------------------------------------
namespace create{

	// Functions
	bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

} // end of namespace create

#endif /*CREATE_H_*/
