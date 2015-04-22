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
#include "demag.hpp"
#include "voronoi.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "random.hpp"
#include "vio.hpp"
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
	std::vector <double> MaterialMuSSIArray(0);
	std::vector <zkval_t> MaterialScalarAnisotropyArray(0);
	std::vector <zkten_t> MaterialTensorAnisotropyArray(0);
   std::vector <double> material_second_order_anisotropy_constant_array(0);
   std::vector <double> material_sixth_order_anisotropy_constant_array(0);
   std::vector <double> material_spherical_harmonic_constants_array(0);
	std::vector <double> MaterialCubicAnisotropyArray(0);

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
	cs::unit_cell_size[0] = 3.0;
	cs::unit_cell_size[1] = 3.0;
	cs::unit_cell_size[2] = 3.0;

	cs::system_dimensions[0] = 100.0;
	cs::system_dimensions[1] = 100.0;
	cs::system_dimensions[2] = 100.0;

	cs::particle_scale   = 50.0;
	cs::particle_spacing = 10.0;
	
	cs::particle_creation_parity=0;
	cs::crystal_structure = "sc";

	// Voronoi Variables
	create_voronoi::voronoi_sd=0.1;
	create_voronoi::parity=0;
	
	// Setup Hamiltonian Flags
	sim::hamiltonian_simulation_flags[0] = 1;	/// Exchange
	sim::hamiltonian_simulation_flags[1] = 1;	/// Anisotropy
	sim::hamiltonian_simulation_flags[2] = 1;	/// Applied
	sim::hamiltonian_simulation_flags[3] = 1;	/// Thermal
	sim::hamiltonian_simulation_flags[4] = 0;	/// Dipolar
	
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
	material[0].Jij_matrix_SI[0]=-11.2e-21;
	material[0].mu_s_SI=1.5*9.27400915e-24;
	material[0].Ku1_SI=-4.644e-24;
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
	cs::unit_cell_size[0] = 3.0;
	cs::unit_cell_size[1] = 3.0;
	cs::unit_cell_size[2] = 3.0;

	cs::system_dimensions[0] = 2.0;
	cs::system_dimensions[1] = 2.0;
	cs::system_dimensions[2] = 2.0;

	cs::particle_scale   = 50.0;
	cs::particle_spacing = 10.0;
	
	cs::particle_creation_parity=0;
	cs::crystal_structure = "sc";
	
	// Turn off multi-spin Flags
	sim::hamiltonian_simulation_flags[0] = 0;	/// Exchange
	sim::hamiltonian_simulation_flags[4] = 0;	/// Dipolar

	// MPI Mode (Homogeneous execution)
	//vmpi::mpi_mode=0;
	//mpi_create_variables::mpi_interaction_range=2; // Unit cells
	//mpi_create_variables::mpi_comms_identify=false;

	return EXIT_SUCCESS;
}


// Simple function to check for valid input for hysteresis loop parameters
void check_hysteresis_loop_parameters(){

   // Only applies to hysteresis loop programs, all others return
   if(sim::program!=12) return;

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
	
	// Ensure H vector is unit length
	// **RE edit 21.11.12 - no longer necessary as value checked on user input**
	//double mod_H=1.0/sqrt(sim::H_vec[0]*sim::H_vec[0]+sim::H_vec[1]*sim::H_vec[1]+sim::H_vec[2]*sim::H_vec[2]);
	//sim::H_vec[0]*=mod_H;
	//sim::H_vec[1]*=mod_H;
	//sim::H_vec[2]*=mod_H;

	// Calculate moment, magnetisation, and anisotropy constants
	/*for(int mat=0;mat<mp::num_materials;mat++){
		double V=cs::unit_cell_size[0]*cs::unit_cell_size[1]*cs::unit_cell_size[2];
		// Set magnetisation from mu_s and a
		if(material[mat].moment_flag==true){
			//material[mat].magnetisation=num_atoms_per_unit_cell*material[mat].mu_s_SI/V;
		}
		// Set mu_s from magnetisation and a
		else {
			//material[mat].mu_s_SI=material[mat].magnetisation*V/num_atoms_per_unit_cell;
		}
		// Set K as energy/atom
		if(material[mat].anis_flag==false){
			material[mat].Ku1_SI=material[mat].Ku1_SI*V/num_atoms_per_unit_cell;
			std::cout << "setting " << material[mat].Ku1_SI << std::endl;
		}
	}*/
	const string blank="";

   // Check for symmetry of exchange matrix
   for(int mi = 0; mi < mp::num_materials; mi++){

      for(int mj = 0; mj < mp::num_materials; mj++){

         // Check for non-zero value (avoids divide by zero)
         if(fabs(material[mi].Jij_matrix_SI[mj]) > 0.0){

            // Calculate ratio of i->j / j-> exchange constants
            double ratio = material[mj].Jij_matrix_SI[mi]/material[mi].Jij_matrix_SI[mj];

            // Check that ratio ~ 1.0 for symmetric exchange interactions
            if( (ratio < 0.99999) || (ratio > 1.00001) ){

               // Error found - report to user and terminate program
               terminaltextcolor(RED);
                  std::cerr << "Error! Non-symmetric exchange interactions for materials " << mi+1 << " and " << mj+1 << ". Exiting" << std::endl;
               terminaltextcolor(WHITE);

               zlog << zTs() << "Error! Non-symmetric exchange interactions for materials " << mi+1 << " and " << mj+1 << std::endl;
               zlog << zTs() << "\tmaterial[" << mi+1 << "]:exchange-matrix[" << mj+1 << "] = " << material[mi].Jij_matrix_SI[mj] << std::endl;
               zlog << zTs() << "\tmaterial[" << mj+1 << "]:exchange-matrix[" << mi+1 << "] = " << material[mj].Jij_matrix_SI[mi] << std::endl;
               zlog << zTs() << "\tThe definition of Heisenberg exchange requires that these values are the same. Exiting." << std::endl;

               err::vexit();

            }
         }
      }
   }

	// Set derived material parameters
	for(int mat=0;mat<mp::num_materials;mat++){
		mp::material[mat].one_oneplusalpha_sq   = -mp::material[mat].gamma_rel/(1.0+mp::material[mat].alpha*mp::material[mat].alpha);
		mp::material[mat].alpha_oneplusalpha_sq =  mp::material[mat].alpha*mp::material[mat].one_oneplusalpha_sq;
			
		for(int j=0;j<mp::num_materials;j++){
			material[mat].Jij_matrix[j]				= mp::material[mat].Jij_matrix_SI[j]/mp::material[mat].mu_s_SI;
		}
		mp::material[mat].Ku									= mp::material[mat].Ku1_SI/mp::material[mat].mu_s_SI;
      mp::material[mat].Ku2                        = mp::material[mat].Ku2_SI/mp::material[mat].mu_s_SI;
      mp::material[mat].Ku3                        = mp::material[mat].Ku3_SI/mp::material[mat].mu_s_SI;
      mp::material[mat].Klatt                      = mp::material[mat].Klatt_SI/mp::material[mat].mu_s_SI;
		mp::material[mat].Kc									= mp::material[mat].Kc1_SI/mp::material[mat].mu_s_SI;
		mp::material[mat].Ks									= mp::material[mat].Ks_SI/mp::material[mat].mu_s_SI;
		mp::material[mat].H_th_sigma						= sqrt(2.0*mp::material[mat].alpha*1.3806503e-23/
																  (mp::material[mat].mu_s_SI*mp::material[mat].gamma_rel*dt));

      // Rename un-named materials with material id
      std::string defname="material#n";
      if(mp::material[mat].name==defname){
         std::stringstream newname;
         newname << "material" << mat+1;
         mp::material[mat].name=newname.str();
      }

      // initialise lattice anisotropy initialisation
      if(sim::lattice_anisotropy_flag==true) mp::material[mat].lattice_anisotropy.set_interpolation_table();

      // output interpolated data to file
      //mp::material[mat].lattice_anisotropy.output_interpolated_function(mat);

	}
		// Check for which anisotropy function(s) are to be used		
		if(sim::TensorAnisotropy==true){
			sim::UniaxialScalarAnisotropy=false; // turn off scalar anisotropy calculation
			// loop over materials and convert all scalar anisotropy to tensor (along z)
			for(int mat=0;mat<mp::num_materials; mat++){
				
				const double one_o_mu=1.0/mp::material[mat].mu_s_SI;

				// If tensor is unset
				if(mp::material.at(mat).KuVec_SI.size()==0){
					const double ex = mp::material.at(mat).UniaxialAnisotropyUnitVector.at(0);
					const double ey = mp::material.at(mat).UniaxialAnisotropyUnitVector.at(1);
					const double ez = mp::material.at(mat).UniaxialAnisotropyUnitVector.at(2);
					mp::material.at(mat).KuVec.push_back(mp::material[mat].Ku*ex*ex);
					mp::material.at(mat).KuVec.push_back(mp::material[mat].Ku*ex*ey);
					mp::material.at(mat).KuVec.push_back(mp::material[mat].Ku*ex*ez);

					mp::material.at(mat).KuVec.push_back(mp::material[mat].Ku*ey*ex);
					mp::material.at(mat).KuVec.push_back(mp::material[mat].Ku*ey*ey);
					mp::material.at(mat).KuVec.push_back(mp::material[mat].Ku*ey*ez);

					mp::material.at(mat).KuVec.push_back(mp::material[mat].Ku*ez*ex);
					mp::material.at(mat).KuVec.push_back(mp::material[mat].Ku*ez*ey);
					mp::material.at(mat).KuVec.push_back(mp::material[mat].Ku*ez*ez);
				}
				else if(mp::material.at(mat).KuVec_SI.size()==9){
					mp::material.at(mat).KuVec.push_back(mp::material.at(mat).KuVec_SI.at(0)*one_o_mu);
					mp::material.at(mat).KuVec.push_back(mp::material.at(mat).KuVec_SI.at(1)*one_o_mu);
					mp::material.at(mat).KuVec.push_back(mp::material.at(mat).KuVec_SI.at(2)*one_o_mu);

					mp::material.at(mat).KuVec.push_back(mp::material.at(mat).KuVec_SI.at(3)*one_o_mu);
					mp::material.at(mat).KuVec.push_back(mp::material.at(mat).KuVec_SI.at(4)*one_o_mu);
					mp::material.at(mat).KuVec.push_back(mp::material.at(mat).KuVec_SI.at(5)*one_o_mu);

					mp::material.at(mat).KuVec.push_back(mp::material.at(mat).KuVec_SI.at(6)*one_o_mu);
					mp::material.at(mat).KuVec.push_back(mp::material.at(mat).KuVec_SI.at(7)*one_o_mu);
					mp::material.at(mat).KuVec.push_back(mp::material.at(mat).KuVec_SI.at(8)*one_o_mu);
				}
			}
		}
		
		// Unroll anisotropy values for speed
		if(sim::UniaxialScalarAnisotropy==true){
			zlog << zTs() << "Setting scalar uniaxial anisotropy." << std::endl;
			// Set global anisotropy type
			sim::AnisotropyType=0;
			MaterialScalarAnisotropyArray.resize(mp::num_materials);
			for(int mat=0;mat<mp::num_materials; mat++) MaterialScalarAnisotropyArray[mat].K=mp::material[mat].Ku;
		}
		else if(sim::TensorAnisotropy==true){
			zlog << zTs() << "Setting tensor uniaxial anisotropy." << std::endl;
			// Set global anisotropy type
			sim::AnisotropyType=1;
			MaterialTensorAnisotropyArray.resize(mp::num_materials);
			for(int mat=0;mat<mp::num_materials; mat++){
				MaterialTensorAnisotropyArray[mat].K[0][0]=mp::material.at(mat).KuVec.at(0);
				MaterialTensorAnisotropyArray[mat].K[0][1]=mp::material.at(mat).KuVec.at(1);
				MaterialTensorAnisotropyArray[mat].K[0][2]=mp::material.at(mat).KuVec.at(2);

				MaterialTensorAnisotropyArray[mat].K[1][0]=mp::material.at(mat).KuVec.at(3);
				MaterialTensorAnisotropyArray[mat].K[1][1]=mp::material.at(mat).KuVec.at(4);
				MaterialTensorAnisotropyArray[mat].K[1][2]=mp::material.at(mat).KuVec.at(5);

				MaterialTensorAnisotropyArray[mat].K[2][0]=mp::material.at(mat).KuVec.at(6);
				MaterialTensorAnisotropyArray[mat].K[2][1]=mp::material.at(mat).KuVec.at(7);
				MaterialTensorAnisotropyArray[mat].K[2][2]=mp::material.at(mat).KuVec.at(8);

			}
		}
      // Unroll second order uniaxial anisotropy values for speed
      if(sim::second_order_uniaxial_anisotropy==true){
         zlog << zTs() << "Setting scalar second order uniaxial anisotropy." << std::endl;
         mp::material_second_order_anisotropy_constant_array.resize(mp::num_materials);
         for(int mat=0;mat<mp::num_materials; mat++) mp::material_second_order_anisotropy_constant_array.at(mat)=mp::material[mat].Ku2;
      }
	  // Unroll sixth order uniaxial anisotropy values for speed
      if(sim::second_order_uniaxial_anisotropy==true){
         zlog << zTs() << "Setting scalar sixth order uniaxial anisotropy." << std::endl;
         mp::material_sixth_order_anisotropy_constant_array.resize(mp::num_materials);
         for(int mat=0;mat<mp::num_materials; mat++) mp::material_sixth_order_anisotropy_constant_array.at(mat)=mp::material[mat].Ku3;
      }
      // Unroll spherical harmonic anisotropy constants for speed
      if(sim::spherical_harmonics==true){
         zlog << zTs() << "Setting spherical harmonics for uniaxial anisotropy" << std::endl;
         mp::material_spherical_harmonic_constants_array.resize(3*mp::num_materials);
         for(int mat=0; mat<mp::num_materials; mat++){
            mp::material_spherical_harmonic_constants_array.at(3*mat+0)=mp::material[mat].sh2/mp::material[mat].mu_s_SI;
            mp::material_spherical_harmonic_constants_array.at(3*mat+1)=mp::material[mat].sh4/mp::material[mat].mu_s_SI;
            mp::material_spherical_harmonic_constants_array.at(3*mat+2)=mp::material[mat].sh6/mp::material[mat].mu_s_SI;
         }
      }
      // Unroll cubic anisotropy values for speed
		if(sim::CubicScalarAnisotropy==true){
			zlog << zTs() << "Setting scalar cubic anisotropy." << std::endl;
			MaterialCubicAnisotropyArray.resize(mp::num_materials);
			for(int mat=0;mat<mp::num_materials; mat++) MaterialCubicAnisotropyArray.at(mat)=mp::material[mat].Kc;
		}

		// Loop over materials to check for invalid input and warn appropriately
		for(int mat=0;mat<mp::num_materials;mat++){
			const double lmin=material[mat].min;
			const double lmax=material[mat].max;
			for(int nmat=0;nmat<mp::num_materials;nmat++){
				if(nmat!=mat){
					double min=material[nmat].min;
					double max=material[nmat].max;
					if(((lmin>min) && (lmin<max)) || ((lmax>min) && (lmax<max))){
						terminaltextcolor(RED);
						std::cerr << "Warning: Overlapping material heights found. Check log for details." << std::endl;
						terminaltextcolor(WHITE);
						zlog << zTs() << "Warning: material " << mat+1 << " overlaps material " << nmat+1 << "." << std::endl;
						zlog << zTs() << "If you have defined geometry then this may be OK, or possibly you meant to specify alloy keyword instead." << std::endl;
						zlog << zTs() << "----------------------------------------------------" << std::endl;
						zlog << zTs() << "  Material "<< mat+1 << ":minimum-height = " << lmin << std::endl;
						zlog << zTs() << "  Material "<< mat+1 << ":maximum-height = " << lmax << std::endl;
						zlog << zTs() << "  Material "<< nmat+1 << ":minimum-height = " << min << std::endl;
						zlog << zTs() << "  Material "<< nmat+1 << ":maximum-height = " << max << std::endl;
					}
				}
			}
		}
	
	return EXIT_SUCCESS;
}

} // end of namespace mp
