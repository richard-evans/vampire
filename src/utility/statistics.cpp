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
///
/// @file
/// @brief Contains functions and variables for calcuation of statistics.
///
/// @details Includes the following functions:
/// \li mag_m
///
/// @section notes Implementation Notes
/// None
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    11/01/2010
/// @internal
///	Created:		11/01/2010
///	Revision:	  ---
///=====================================================================================
///
// Headers
#include "anisotropy.hpp"
#include "atoms.hpp"
#include "exchange.hpp"
#include "gpu.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "sim.hpp"
#include "stats.hpp"

#include <cmath>
#include <iostream>

//Function prototypes
int calculate_spin_fields(const int,const int);
int calculate_external_fields(const int,const int);

/// @namespace stats
/// @brief Variables and functions for calculation of system statistics.
///
/// @internal
///=====================================================================================
///
namespace stats
{
   int num_atoms;				/// Number of atoms for statistic purposes
   double inv_num_atoms;	/// 1.0/num_atoms
   double max_moment;		/// Total Maximum moment

	// torque calculation
	bool calculate_torque=false;
	double total_system_torque[3]={0.0,0.0,0.0};
	double total_mean_system_torque[3]={0.0,0.0,0.0};

	std::vector <double> sublattice_mean_torque_x_array(0);
	std::vector <double> sublattice_mean_torque_y_array(0);
	std::vector <double> sublattice_mean_torque_z_array(0);

	double torque_data_counter=0.0;

	// function prototypes
	void system_torque();
	//void system_energy();

	bool is_initialised=false;

	double data_counter=0.0;		/// number of data points for averaging

// Namespace Functions
/// @brief Calculates total and sublattice, normalised and actual, magnetic moments of the system.
///
/// @details For single materials an optimised version of this routine is used.
///
/// @section notes Implementation Notes
/// None.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author		Richard Evans, rfle500@york.ac.uk
/// @version	1.0
/// @date		11/01/2010
/// @todo		Implement for more than one material - Urgent!
///
/// @return exit_success
///
/// @internal
///	Created:		11/01/2010
///	Revision:	  ---
///=====================================================================================
///
int mag_m(){

   //------------------------------------------------------------------
   // Calculate number and inverse number of moments for normalisation
   //------------------------------------------------------------------
   #ifdef MPICF
      stats::num_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
   #else
      stats::num_atoms = atoms::num_atoms;
   #endif

   stats::inv_num_atoms = 1.0/double(stats::num_atoms);

   if(!stats::is_initialised){

      // Calculate number of moments in each sublattice
      for(int atom=0;atom<stats::num_atoms;atom++){
         int mat=atoms::type_array[atom];
         // add to max_moment
         stats::max_moment+=mp::material[mat].mu_s_SI;
      }

      // Calculate global moment for all CPUs
      #ifdef MPICF
         MPI_Allreduce(MPI_IN_PLACE,&stats::max_moment,1,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
      #endif

      // Resize arrays
      stats::sublattice_mean_torque_x_array.resize(mp::num_materials,0.0);
      stats::sublattice_mean_torque_y_array.resize(mp::num_materials,0.0);
      stats::sublattice_mean_torque_z_array.resize(mp::num_materials,0.0);

      // Set initilaised flag to true
      stats::is_initialised=true;

   }

   // update statistics - need to eventually replace mag_m() with stats::update()...
   stats::update(atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array, atoms::m_spin_array, atoms::type_array, sim::temperature);

   // optionally calculate system torque
   if(stats::calculate_torque==true) stats::system_torque();

   // increment data counter
   stats::data_counter+=1.0;

   return EXIT_SUCCESS;
}

/// @brief Resets mean magnetisation and counter.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author		Richard Evans, richard.evans@york.ac.uk
/// @version	1.1
/// @date		14/09/2011
///
/// @internal
///	Created:		11/01/2010
///	Revision:	  ---
///=====================================================================================
///
void mag_m_reset(){
	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "stats::mag_m_reset() has been called" << std::endl;}

	// if stats not initialised then call mag_m() to do so.
	if(!stats::is_initialised) stats::mag_m();

   // reset statistics - need to eventually replace mag_m_reset() with stats::reset()...
   stats::reset();

	stats::data_counter=0.0;

	stats::total_mean_system_torque[0]=0.0;
	stats::total_mean_system_torque[1]=0.0;
	stats::total_mean_system_torque[2]=0.0;

	for(int mat=0;mat<mp::num_materials;mat++){
		stats::sublattice_mean_torque_x_array[mat]=0.0;
		stats::sublattice_mean_torque_y_array[mat]=0.0;
		stats::sublattice_mean_torque_z_array[mat]=0.0;
	}
	stats::torque_data_counter=0.0;

}

double max_torque(){
  ///================================================================================================
  ///
 ///                                 subroutine torque
  ///
  ///                      Calculates total torque on the system
  ///
  ///================================================================================================

	double max_torque=0.0;
	double mag_torque;
	double torque[3];

	//------------------------------------------------
	// Recalculate net fields
	//------------------------------------------------

	calculate_spin_fields(0,num_atoms);
	calculate_external_fields(0,num_atoms);

	for(int atom=0;atom<num_atoms;atom++){

		// Store local spin in Sand local field in H
		const double S[3] = {atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};
		const double H[3] = {atoms::x_total_spin_field_array[atom]+atoms::x_total_external_field_array[atom],
									atoms::y_total_spin_field_array[atom]+atoms::y_total_external_field_array[atom],
									atoms::z_total_spin_field_array[atom]+atoms::z_total_external_field_array[atom]};

		torque[0] = S[1]*H[2]-S[2]*H[1];
		torque[1] = S[2]*H[0]-S[0]*H[2];
		torque[2] = S[0]*H[1]-S[1]*H[0];

		mag_torque = sqrt(torque[0]*torque[0] + torque[1]*torque[1] + torque[2]*torque[2]);

		if(mag_torque>max_torque){
			max_torque = mag_torque;
		}

	}

   // find max torque on all nodes
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE,&max_torque,1,MPI_DOUBLE,MPI_MAX, MPI_COMM_WORLD);
   #endif

  return max_torque;

}

///
/// @brief      Calculates the torque on the system
///
/// Calculates the instantaneous value of the system torque T = sum (Si x Hi)
///
/// @return     void
///
void system_torque(){

	// parallel version, multiple materials, initialisation

	double torque[3]={0.0,0.0,0.0};

	// calculate net fields
	calculate_spin_fields(0,atoms::num_atoms);
	calculate_external_fields(0,atoms::num_atoms);

	for(int atom=0;atom<stats::num_atoms;atom++){

		// get atomic moment
		const int imat=atoms::type_array[atom];
		const double mu = mp::material[imat].mu_s_SI;

		// Store local spin in Sand local field in H
		const double S[3] = {atoms::x_spin_array[atom]*mu,atoms::y_spin_array[atom]*mu,atoms::z_spin_array[atom]*mu};
		const double H[3] = {atoms::x_total_spin_field_array[atom]+atoms::x_total_external_field_array[atom],
									atoms::y_total_spin_field_array[atom]+atoms::y_total_external_field_array[atom],
									atoms::z_total_spin_field_array[atom]+atoms::z_total_external_field_array[atom]};

		torque[0] += S[1]*H[2]-S[2]*H[1];
		torque[1] += S[2]*H[0]-S[0]*H[2];
		torque[2] += S[0]*H[1]-S[1]*H[0];

		stats::sublattice_mean_torque_x_array[imat]+=S[1]*H[2]-S[2]*H[1];
		stats::sublattice_mean_torque_y_array[imat]+=S[2]*H[0]-S[0]*H[2];
		stats::sublattice_mean_torque_z_array[imat]+=S[0]*H[1]-S[1]*H[0];
	}

	// reduce torque on all nodes
	#ifdef MPICF
		MPI_Allreduce(MPI_IN_PLACE,&torque[0],3,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE,&stats::sublattice_mean_torque_x_array[0],mp::num_materials,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE,&stats::sublattice_mean_torque_y_array[0],mp::num_materials,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE,&stats::sublattice_mean_torque_z_array[0],mp::num_materials,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
	#endif

	// Set stats values
	stats::total_system_torque[0]=torque[0];
	stats::total_system_torque[1]=torque[1];
	stats::total_system_torque[2]=torque[2];

	stats::total_mean_system_torque[0]+=torque[0];
	stats::total_mean_system_torque[1]+=torque[1];
	stats::total_mean_system_torque[2]+=torque[2];

	stats::torque_data_counter+=1.0;
	return;

}

} // End of Namespace
