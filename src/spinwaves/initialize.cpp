//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sergiu Ruta 2022. All rights reserved.
//
//   Email: sergiu.ruta@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "spinwaves.hpp"

// sw module headers
#include "internal.hpp"
#include "iostream"
#include <iomanip>
#include <string>
#include <sstream>
#include "vmpi.hpp"

//sergiu for SW
#include "unitcell.hpp"
#include "errors.hpp"
#include "sim.hpp"
#include "vio.hpp"
#include "vector"
#include <cmath>
#include "fstream"
#include "atoms.hpp"
#include "program.hpp"



namespace spinwaves{
	//   std::vector <double> Skx_FFT_array;
	//   std::vector <double> Sky_FFT_array;
	//   std::vector <double> Skz_FFT_array;

	//----------------------------------------------------------------------------
	// Function to initialize sw module
	//----------------------------------------------------------------------------
	void initialize(const double system_dimensions_x,
					const double system_dimensions_y,
					const double system_dimensions_z,
					const double total_num_unit_cells_x,
					const double total_num_unit_cells_y,
					const double total_num_unit_cells_z,
					const double unit_cell_size_x,
					const double unit_cell_size_y,
					const double unit_cell_size_z,
					const std::vector <unitcell::atom_t>& atom,
					const std::vector<double>& atom_coords_x,
					const std::vector<double>& atom_coords_y,
					const std::vector<double>& atom_coords_z){

		//-------------------------------------------------------------------------------------
		// Check if spinwave calculation enabled, if not do nothing
		//-------------------------------------------------------------------------------------
		if(program::program!=74) return;

		// check spectrum values in input file are in agreement with the spinwaves:number-of-spectrums
		spinwaves::internal::check_numbering_of_spectrums();

		int Na=atom.size();
		const double toll_Fikj=1e-5; // tollerance for the structure factor values smaller than tollerance are considered 0

		// number of steps for sw calculations
		int numsteps = (sim::total_time / sim::partial_time);

		// if caclulating spinwaves for a specific material, we need to apply a mask
		spinwaves::internal::calculate_material_mask();

		// JRH determine whether to use path for a predefined crystal or whether path has been specified by user.
		spinwaves::internal::determine_path();

		// JRH calculate the values of $k$ at which to calculate the DSF for the user specified system dimensions.
		if (internal::filetype == "specific-k") spinwaves::internal::determine_kpoints_from_user_specific_k(system_dimensions_x, system_dimensions_y, system_dimensions_z, unit_cell_size_x, unit_cell_size_y, unit_cell_size_z);
		if (internal::filetype == "path") spinwaves::internal::determine_kpoints_from_user_high_sym_path(system_dimensions_x, system_dimensions_y, system_dimensions_z, unit_cell_size_x, unit_cell_size_y, unit_cell_size_z);

		// initialise arrays based on kpoints calculated above
		spinwaves::internal::initialise_arrays();

		// determine prefactor that will be used in fourier transform. sin(k_x*r_x) etc.
		spinwaves::internal::calculate_fourier_prefactor(atom_coords_x, atom_coords_y, atom_coords_z);

		// determine which component of spin to calculate spinwave dispersion from
		spinwaves::internal::determine_spin_component();

      	// output a file containing frequencies
    	spinwaves::internal::save_frequencies();

		std::cout<< "SW initialisation completed."<<std::endl;
		std::cout<<" "<<std::endl;


		return;

	}

} // end of sw namespace
//
