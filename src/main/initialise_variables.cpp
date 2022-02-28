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

// Headers
#include "errors.hpp"
#include "dipole.hpp"
#include "voronoi.hpp"
#include "material.hpp"
#include "program.hpp"
#include "sim.hpp"
#include "random.hpp"
#include "vio.hpp"
#include "unitcell.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

#include <cmath>
#include <iostream>
#include <sstream>
//==========================================================
// Namespace material_parameters
//==========================================================
namespace mp{
	//----------------------------------
	// Material Container
	//----------------------------------

	//const int max_materials=100;

	int num_materials=1;


	std::vector <materials_t> material(1);



	//----------------------------------
	//Input Integration parameters
	//----------------------------------
	double dt_SI;
	double gamma_SI = 1.76E11;

	//----------------------------------
	//Derived Integration parameters
	//----------------------------------
	double dt;
	double half_dt;

	// Unrolled material parameters for speed
	std::vector <double> mu_s_array;

///
/// @brief Function to initialise program variables prior to system creation.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    19/01/2010
///
/// @param[in] infile Main input file name for system initialisation
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		19/01/2010
///	Revision:	  ---
///=====================================================================================
///
int initialise(std::string const infile){
	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "initialise_variables has been called" << std::endl;}

	if(vmpi::my_rank==0){
		std::cout << "================================================================================" << std::endl;
		std::cout << "Initialising system variables" << std::endl;
	}

	// Setup default system settings
	mp::default_system();

	// Read values from input files
	int iostat = vin::read(infile);
	if(iostat==EXIT_FAILURE){
		terminaltextcolor(RED);
		std::cerr << "Error - input file \'" << infile << "\' not found, exiting" << std::endl;
		terminaltextcolor(WHITE);
		err::vexit();
	}

	// Print out material properties
	//mp::material[0].print();

	// Check for keyword parameter overide
	if(cs::single_spin==true){
		mp::single_spin_system();
	}

	// Set derived system parameters
	mp::set_derived_parameters();

	// Return
	return EXIT_SUCCESS;
}

int default_system(){

	// Initialise system creation flags to zero
	for (int i=0;i<10;i++){
		cs::system_creation_flags[i] = 0;
		sim::hamiltonian_simulation_flags[i] = 0;
	}

	// Set system dimensions !Angstroms

	cs::system_dimensions[0] = 100.0;
	cs::system_dimensions[1] = 100.0;
	cs::system_dimensions[2] = 100.0;

	cs::particle_scale   = 50.0;
	cs::particle_spacing = 10.0;

   cs::particle_creation_parity=0;
   uc::set_crystal_structure_to_simple_cubic();

	// Voronoi Variables
	create_voronoi::voronoi_sd=0.1;
	create_voronoi::parity=0;

	// Setup Hamiltonian Flags
	sim::hamiltonian_simulation_flags[0] = 1;	/// Exchange
	sim::hamiltonian_simulation_flags[1] = 1;	/// Anisotropy
	sim::hamiltonian_simulation_flags[2] = 1;	/// Applied
	sim::hamiltonian_simulation_flags[3] = 1;	/// Thermal

	//Integration parameters
	dt_SI = 1.0e-15;	// seconds
	dt = dt_SI*mp::gamma_SI; // Must be set before Hth
	half_dt = 0.5*dt;

	//------------------------------------------------------------------------------
	// Material Definitions
	//------------------------------------------------------------------------------
	num_materials=1;
	material.resize(num_materials);

	//-------------------------------------------------------
	// Material 0
	//-------------------------------------------------------
	material[0].name="Co";
	material[0].alpha=0.1;
	material[0].mu_s_SI=1.5*9.27400915e-24;
	material[0].gamma_rel=1.0;
	material[0].element="Ag ";

	// Disable Error Checking
	err::check=false;

	// Initialise random number generator
	mtrandom::grnd.seed(2106975519);

	return EXIT_SUCCESS;
}

int single_spin_system(){

	// Reset system creation flags to zero
	for (int i=0;i<10;i++){
		cs::system_creation_flags[i] = 0;
	}

	// Set system dimensions !Angstroms

	cs::system_dimensions[0] = 2.0;
	cs::system_dimensions[1] = 2.0;
	cs::system_dimensions[2] = 2.0;

	cs::particle_scale   = 50.0;
	cs::particle_spacing = 10.0;

	cs::particle_creation_parity=0;
	uc::set_crystal_structure_to_simple_cubic();

	// Turn off multi-spin Flags
	sim::hamiltonian_simulation_flags[0] = 0;	/// Exchange

	// MPI Mode (Homogeneous execution)
	//vmpi::mpi_mode=0;
	//mpi_create_variables::mpi_interaction_range=2; // Unit cells
	//mpi_create_variables::mpi_comms_identify=false;

	return EXIT_SUCCESS;
}


// Simple function to check for valid input for hysteresis loop parameters
void check_hysteresis_loop_parameters(){

   // Only applies to hysteresis loop programs, all others return
   if(program::program!=12) return;

   double min=sim::Hmin;
   double max=sim::Hmax;
   double inc=sim::Hinc;

   // + + +
   if(min>=0 && max>=0 && inc>0){
      if(max<min){
         if(vmpi::my_rank==0){
			terminaltextcolor(RED);
            std::cout << "Error in hysteresis-loop parameters:" << std::endl;
            std::cout << "\t sim:minimum-applied-field-strength = " << min << std::endl;
            std::cout << "\t sim:maximum-applied-field-strength = " << max << std::endl;
            std::cout << "\t sim:applied-field-strength-increment = " << inc << std::endl;
            std::cout << "Minimum and maximum fields are both positive, but minimum > maximum with a positive increment, causing an infinite loop. Exiting." << std::endl;
			terminaltextcolor(WHITE);
            zlog << zTs() << "Error in hysteresis-loop parameters:" << std::endl;
            zlog << zTs() << "\t sim:minimum-applied-field-strength = " << min << std::endl;
            zlog << zTs() << "\t sim:maximum-applied-field-strength = " << max << std::endl;
            zlog << zTs() << "\t sim:applied-field-strength-increment = " << inc << std::endl;
            zlog << zTs() << "Minimum and maximum fields are both positive, but minimum > maximum with a positive increment, causing an infinite loop. Exiting." << std::endl;
            err::vexit();
         }
      }
   }
   // + + -
   else if(min>=0 && max>=0 && inc<0){
      if(max>min){
         if(vmpi::my_rank==0){
			terminaltextcolor(RED);
            std::cout << "Error in hysteresis-loop parameters:" << std::endl;
            std::cout << "\t sim:minimum-applied-field-strength = " << min << std::endl;
            std::cout << "\t sim:maximum-applied-field-strength = " << max << std::endl;
            std::cout << "\t sim:applied-field-strength-increment = " << inc << std::endl;
            std::cout << "Minimum and maximum fields are both positive, but maximum > minimum with a negative increment, causing an infinite loop. Exiting." << std::endl;
            terminaltextcolor(WHITE);
			zlog << zTs() << "Error in hysteresis-loop parameters:" << std::endl;
            zlog << zTs() << "\t sim:minimum-applied-field-strength = " << min << std::endl;
            zlog << zTs() << "\t sim:maximum-applied-field-strength = " << max << std::endl;
            zlog << zTs() << "\t sim:applied-field-strength-increment = " << inc << std::endl;
            zlog << zTs() << "Minimum and maximum fields are both positive, but maximum > minimum with a negative increment, causing an infinite loop. Exiting." << std::endl;
            err::vexit();
         }
      }
   }
   // + - +
   else if(min>=0 && max<0 && inc>0){
      if(vmpi::my_rank==0){
		 terminaltextcolor(RED);
         std::cout << "Error in hysteresis-loop parameters:" << std::endl;
         std::cout << "\t sim:minimum-applied-field-strength = " << min << std::endl;
         std::cout << "\t sim:maximum-applied-field-strength = " << max << std::endl;
         std::cout << "\t sim:applied-field-strength-increment = " << inc << std::endl;
         std::cout << "Minimum field is positive and maximum field is negative with a positive increment, causing an infinite loop. Exiting." << std::endl;
         terminaltextcolor(WHITE);
		 zlog << zTs() << "Error in hysteresis-loop parameters:" << std::endl;
         zlog << zTs() << "\t sim:minimum-applied-field-strength = " << min << std::endl;
         zlog << zTs() << "\t sim:maximum-applied-field-strength = " << max << std::endl;
         zlog << zTs() << "\t sim:applied-field-strength-increment = " << inc << std::endl;
         zlog << zTs() << "Minimum field is positive and maximum field is negative with a positive increment, causing an infinite loop. Exiting." << std::endl;
         err::vexit();
      }
   }
   // - + -
   else if(min<0 && max>=0 && inc<0){
      if(vmpi::my_rank==0){
		 terminaltextcolor(RED);
         std::cout << "Error in hysteresis-loop parameters:" << std::endl;
         std::cout << "\t sim:minimum-applied-field-strength = " << min << std::endl;
         std::cout << "\t sim:maximum-applied-field-strength = " << max << std::endl;
         std::cout << "\t sim:applied-field-strength-increment = " << inc << std::endl;
         std::cout << "Minimum field is negative and maximum field is positive with a negative increment, causing an infinite loop. Exiting." << std::endl;
         terminaltextcolor(WHITE);
		 zlog << zTs() << "Error in hysteresis-loop parameters:" << std::endl;
         zlog << zTs() << "\t sim:minimum-applied-field-strength = " << min << std::endl;
         zlog << zTs() << "\t sim:maximum-applied-field-strength = " << max << std::endl;
         zlog << zTs() << "\t sim:applied-field-strength-increment = " << inc << std::endl;
         zlog << zTs() << "Minimum field is negative and maximum field is positive with a negative increment, causing an infinite loop. Exiting." << std::endl;
         err::vexit();
      }
   }
   // - - -
   else if(min<0 && max<0 && inc<0){
      if(max>min){
         if(vmpi::my_rank==0){
			terminaltextcolor(RED);
            std::cout << "Error in hysteresis-loop parameters:" << std::endl;
            std::cout << "\t sim:minimum-applied-field-strength = " << min << std::endl;
            std::cout << "\t sim:maximum-applied-field-strength = " << max << std::endl;
            std::cout << "\t sim:applied-field-strength-increment = " << inc << std::endl;
            std::cout << "Minimum and maximum fields are both negative, but minimum < maximum with a negative increment, causing an infinite loop. Exiting." << std::endl;
            terminaltextcolor(WHITE);
			zlog << zTs() << "Error in hysteresis-loop parameters:" << std::endl;
            zlog << zTs() << "\t sim:minimum-applied-field-strength = " << min << std::endl;
            zlog << zTs() << "\t sim:maximum-applied-field-strength = " << max << std::endl;
            zlog << zTs() << "\t sim:applied-field-strength-increment = " << inc << std::endl;
            zlog << zTs() << "Minimum and maximum fields are both negative, but minimum < maximum with a negative increment, causing an infinite loop. Exiting." << std::endl;
            err::vexit();
         }
      }
   }
   // - - +
   else if(min<0 && max<0 && inc>0){
      if(max<min){
         if(vmpi::my_rank==0){
			terminaltextcolor(RED);
            std::cout << "Error in hysteresis-loop parameters:" << std::endl;
            std::cout << "\t sim:minimum-applied-field-strength = " << min << std::endl;
            std::cout << "\t sim:maximum-applied-field-strength = " << max << std::endl;
            std::cout << "\t sim:applied-field-strength-increment = " << inc << std::endl;
            std::cout << "Minimum and maximum fields are both negative, but maximum < minimum with a positive increment, causing an infinite loop. Exiting." << std::endl;
            terminaltextcolor(WHITE);
			zlog << zTs() << "Error in hysteresis-loop parameters:" << std::endl;
            zlog << zTs() << "\t sim:minimum-applied-field-strength = " << min << std::endl;
            zlog << zTs() << "\t sim:maximum-applied-field-strength = " << max << std::endl;
            zlog << zTs() << "\t sim:applied-field-strength-increment = " << inc << std::endl;
            zlog << zTs() << "Minimum and maximum fields are both positive, but maximum < minimum with a positive increment, causing an infinite loop. Exiting." << std::endl;
            err::vexit();
         }
      }
   }
   return;

}

int set_derived_parameters(){

	// Set integration constants
	mp::dt = mp::dt_SI*mp::gamma_SI; // Must be set before Hth
	mp::half_dt = 0.5*mp::dt;

	// Check to see if field direction is set by angle
	if(sim::applied_field_set_by_angle){
		sim::H_vec[0]=sin(sim::applied_field_angle_phi*M_PI/180.0)*cos(sim::applied_field_angle_theta*M_PI/180.0);
		sim::H_vec[1]=sin(sim::applied_field_angle_phi*M_PI/180.0)*sin(sim::applied_field_angle_theta*M_PI/180.0);
		sim::H_vec[2]=cos(sim::applied_field_angle_phi*M_PI/180.0);
	}

	// Check for valid particle array offsets
	if(cs::particle_array_offset_x >= cs::system_dimensions[0]){
		terminaltextcolor(RED);
		std::cerr << "Warning: requested particle-array-offset-x is greater than system dimensions." << std::endl;
		std::cerr << "Info: This will probably lead to no particles being created and generate an error." << std::endl;
		terminaltextcolor(WHITE);
		zlog << zTs() << "Warning: requested particle-array-offset-x is greater than system dimensions." << std::endl;
		zlog << zTs() << "Info: This will probably lead to no particles being created and generate an error." << std::endl;
	}
	if(cs::particle_array_offset_y >= cs::system_dimensions[1]){
		terminaltextcolor(RED);
		std::cerr << "Warning: requested particle-array-offset-y is greater than system dimensions." << std::endl;
		std::cerr << "Info: This will probably lead to no particles being created and generate an error." << std::endl;
		terminaltextcolor(WHITE);
		zlog << zTs() << "Warning: requested particle-array-offset-y is greater than system dimensions." << std::endl;
		zlog << zTs() << "Info: This will probably lead to no particles being created and generate an error." << std::endl;
	}

	check_hysteresis_loop_parameters();

	const string blank="";

	// Set derived material parameters
	for(int mat=0;mat<mp::num_materials;mat++){
		mp::material[mat].one_oneplusalpha_sq   = -mp::material[mat].gamma_rel/(1.0+mp::material[mat].alpha*mp::material[mat].alpha);
		mp::material[mat].alpha_oneplusalpha_sq =  mp::material[mat].alpha*mp::material[mat].one_oneplusalpha_sq;
		mp::material[mat].H_th_sigma			    = sqrt(2.0*mp::material[mat].alpha*1.3806503e-23 / (mp::material[mat].mu_s_SI*mp::material[mat].gamma_rel*dt));

      // Rename un-named materials with material id
      std::string defname="material#n";
      if(mp::material[mat].name==defname){
         std::stringstream newname;
         newname << "material" << mat+1;
         mp::material[mat].name=newname.str();
      }

	}

      // Unroll material spin moment values for speed
      mp::mu_s_array.resize(mp::num_materials);
      for(int mat=0;mat<mp::num_materials; mat++) mu_s_array.at(mat)=mp::material[mat].mu_s_SI/9.27400915e-24; // normalise to mu_B

	return EXIT_SUCCESS;
}

} // end of namespace mp
