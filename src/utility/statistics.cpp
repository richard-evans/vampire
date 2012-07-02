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
#include "atoms.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

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
	int num_atoms;				// Number of atoms for statistic purposes
	double inv_num_atoms;	//1.0/num_atoms
	double max_moment;		// Total Maximum moment

	double total_mag_actual[3];	///< Actual magnetisation components
	double total_mag_m_actual;		///< Actual magnitude of total magnetisation
	
	double total_mag_norm[3];	///< Normalised magnetisation components
	double total_mag_m_norm;	///< Normalised magnitude of total magnetisation
	
	double total_mean_mag_m_actual=0.0;
	double total_mean_mag_m_norm=0.0;
	
	std::vector <double> sublattice_mx_array(0);
	std::vector <double> sublattice_my_array(0);
	std::vector <double> sublattice_mz_array(0);
	std::vector <double> sublattice_magm_array(0);
	std::vector <double> sublattice_mean_magm_array(0);
	std::vector <double> sublattice_mom_array(0);
	std::vector <int> sublattice_nm_array(0);
	
	// torque calculation
	bool calculate_torque=false;
	double total_system_torque[3]={0.0,0.0,0.0};
	double total_mean_system_torque[3]={0.0,0.0,0.0};
	
	std::vector <double> sublattice_mean_torque_x_array(0);
	std::vector <double> sublattice_mean_torque_y_array(0);
	std::vector <double> sublattice_mean_torque_z_array(0);
	
	double torque_data_counter=0.0;
	
	// susceptibility (chi) calculation
	bool CalculateSusceptibility=false;
	double MeanChi[3]={0.0,0.0,0.0};
	double MeanChiSquared[3]={0.0,0.0,0.0};
	double MeanChiDataCounter=0.0;
	double ChiAtoms=0;
	
	// function prototypes
	void system_torque();
	void SystemSusceptibility();
	
	bool is_initialised=false;

	double data_counter=0.0;		// number of data points for averaging
	
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
	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "stats::mag_m() has been called" << std::endl;}
	
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
		// Resize arrays
		stats::sublattice_mx_array.resize(mp::num_materials,0.0);
		stats::sublattice_my_array.resize(mp::num_materials,0.0);
		stats::sublattice_mz_array.resize(mp::num_materials,0.0);
		stats::sublattice_magm_array.resize(mp::num_materials,0.0);
		stats::sublattice_mean_magm_array.resize(mp::num_materials,0.0);
		stats::sublattice_mom_array.resize(mp::num_materials,0.0);
		stats::sublattice_nm_array.resize(mp::num_materials,0);
		stats::sublattice_mean_torque_x_array.resize(mp::num_materials,0.0);
		stats::sublattice_mean_torque_y_array.resize(mp::num_materials,0.0);
		stats::sublattice_mean_torque_z_array.resize(mp::num_materials,0.0);
		stats::max_moment=0.0;
		
		// Calculate number of moments in each sublattice
		for(int atom=0;atom<stats::num_atoms;atom++){
			int mat=atoms::type_array[atom];
			
			// count towards moment
			stats::sublattice_nm_array[mat]++;
			// add to max_moment
			stats::max_moment+=mp::material[mat].mu_s_SI;
		}
		
		int nm=0;
		// Calculate size of each moment and check num_moments
		for(int mat=0;mat<mp::num_materials;mat++){
			stats::sublattice_mom_array[mat]=mp::material[mat].mu_s_SI;
			nm+=stats::sublattice_nm_array[mat];
		}
		// Check correct number of moments found
		if(nm != stats::num_atoms){
			std::cerr << "Error in mag_m calculation, missing moments are present:" << stats::num_atoms << " expected, " << nm << " found!" << std::endl;
		}

		// Calculate global moment for all CPUs
		#ifdef MPICF
			MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&stats::max_moment,1,MPI_DOUBLE,MPI_SUM);
			MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&stats::sublattice_nm_array[0],mp::num_materials,MPI_INT,MPI_SUM);
		#endif

		// Output maximum moment to file
		zlog << zTs() << "Maximum spin moment for sample is " << stats::max_moment << " [J/T]" << std::endl;
		zlog << zTs() << "Maximum spin moment per atom is " << (stats::max_moment/double(stats::num_atoms))/9.27400915e-24 << " [mu_B/atom]" << std::endl;
			
		// Set initilaised flag to true
		stats::is_initialised=true;
	}	
	
	// optionally calculate system torque
	if(stats::calculate_torque==true) stats::system_torque();
	
	//optionally calculate system susceptibility
	if(stats::CalculateSusceptibility==true) stats::SystemSusceptibility();
	
	//---------------------------------------------------------
	// For single material, use streamlined version for speed
	//---------------------------------------------------------
	if(material_parameters::num_materials==1){
		double m[3]={0.0,0.0,0.0};

		// Calculate total components
		for(int atom=0;atom<stats::num_atoms;atom++){
			m[0] += atoms::x_spin_array[atom];
			m[1] += atoms::y_spin_array[atom];
			m[2] += atoms::z_spin_array[atom];
		}
		
		#ifdef MPICF
			MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&stats::num_atoms,1,MPI_INT,MPI_SUM);
			MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&m[0],1,MPI_DOUBLE,MPI_SUM);
			MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&m[1],1,MPI_DOUBLE,MPI_SUM);
			MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&m[2],1,MPI_DOUBLE,MPI_SUM);
			//MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&stats::mag_m,1,MPI_DOUBLE,MPI_SUM);
		#endif

		stats::total_mag_norm[0]=m[0]/double(stats::num_atoms);
		stats::total_mag_norm[1]=m[1]/double(stats::num_atoms);
		stats::total_mag_norm[2]=m[2]/double(stats::num_atoms);
		stats::total_mag_m_norm = sqrt(stats::total_mag_norm[0]*stats::total_mag_norm[0] +
												 stats::total_mag_norm[1]*stats::total_mag_norm[1] +
												 stats::total_mag_norm[2]*stats::total_mag_norm[2]);

		// Calculate actual moments
		stats::total_mag_actual[0] = stats::total_mag_norm[0]*stats::max_moment;
		stats::total_mag_actual[1] = stats::total_mag_norm[1]*stats::max_moment;
		stats::total_mag_actual[2] = stats::total_mag_norm[2]*stats::max_moment;
		stats::total_mag_m_actual  = stats::total_mag_m_norm*stats::max_moment;

		// Add to mean magm
		stats::total_mean_mag_m_norm+=stats::total_mag_m_norm;
		stats::total_mean_mag_m_actual+=stats::total_mag_m_actual;
		
		// Copy to sublattice arrays for output
		stats::sublattice_mx_array[0] = stats::total_mag_norm[0];
		stats::sublattice_my_array[0] = stats::total_mag_norm[1];
		stats::sublattice_mz_array[0] = stats::total_mag_norm[2];
		stats::sublattice_magm_array[0] = stats::total_mag_m_norm;
		stats::sublattice_mean_magm_array[0] = stats::total_mean_mag_m_norm;
		
		// Increment data counter		
		stats::data_counter+=1.0;
	}
	else{
	//---------------------------------------------------------
	// For multiple materials, use full version
	//---------------------------------------------------------
		
		// Reinitialise sublattice components
		for(int mat=0;mat<mp::num_materials;mat++){
			stats::sublattice_mx_array[mat]=0.0;
			stats::sublattice_my_array[mat]=0.0;
			stats::sublattice_mz_array[mat]=0.0;
		}

		// Calculate sublattice components
		for(int atom=0;atom<stats::num_atoms;atom++){
			int mat=atoms::type_array[atom];
			stats::sublattice_mx_array[mat]+=atoms::x_spin_array[atom];
			stats::sublattice_my_array[mat]+=atoms::y_spin_array[atom];
			stats::sublattice_mz_array[mat]+=atoms::z_spin_array[atom];
		}

		// Reduce sublattice moments on all CPUs
		#ifdef MPICF
			MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&stats::sublattice_mx_array[0],mp::num_materials,MPI_DOUBLE,MPI_SUM);
			MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&stats::sublattice_my_array[0],mp::num_materials,MPI_DOUBLE,MPI_SUM);
			MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&stats::sublattice_mz_array[0],mp::num_materials,MPI_DOUBLE,MPI_SUM);
		#endif
		
		// initialise total variables
		stats::total_mag_m_actual = 0.0;
		stats::total_mag_actual[0]= 0.0;
		stats::total_mag_actual[1]= 0.0;
		stats::total_mag_actual[2]= 0.0;
		
		// Calculate sublattice magnitude
		for(int mat=0;mat<mp::num_materials;mat++){
			double inm = 1.0/stats::sublattice_nm_array[mat];
			double mx=stats::sublattice_mx_array[mat]*inm;
			double my=stats::sublattice_my_array[mat]*inm;
			double mz=stats::sublattice_mz_array[mat]*inm;

			// Calculate magnitudes
			stats::sublattice_magm_array[mat]=sqrt(mx*mx + my*my + mz*mz);
			stats::sublattice_mean_magm_array[mat]+=stats::sublattice_magm_array[mat];
						
			// Calculate actual total moment
			stats::total_mag_actual[0]+=stats::sublattice_mx_array[mat]*stats::sublattice_mom_array[mat];
			stats::total_mag_actual[1]+=stats::sublattice_my_array[mat]*stats::sublattice_mom_array[mat];
			stats::total_mag_actual[2]+=stats::sublattice_mz_array[mat]*stats::sublattice_mom_array[mat];
			
			// save (un)normalised moments
			stats::sublattice_mx_array[mat]=mx; // /stats::sublattice_magm_array[mat];
			stats::sublattice_my_array[mat]=my; // /stats::sublattice_magm_array[mat];
			stats::sublattice_mz_array[mat]=mz; // /stats::sublattice_magm_array[mat];
		}

		stats::total_mag_norm[0]=stats::total_mag_actual[0]/stats::max_moment;
		stats::total_mag_norm[1]=stats::total_mag_actual[1]/stats::max_moment;
		stats::total_mag_norm[2]=stats::total_mag_actual[2]/stats::max_moment;

		// Calculate total moment
		total_mag_m_norm = sqrt(stats::total_mag_norm[0]*stats::total_mag_norm[0] + 
										stats::total_mag_norm[1]*stats::total_mag_norm[1] +
										stats::total_mag_norm[2]*stats::total_mag_norm[2]);
		
		total_mag_m_actual = sqrt(stats::total_mag_actual[0]*stats::total_mag_actual[0] + 
										  stats::total_mag_actual[1]*stats::total_mag_actual[1] +
										  stats::total_mag_actual[2]*stats::total_mag_actual[2]);

		stats::total_mean_mag_m_norm+=stats::total_mag_m_norm;
		stats::total_mean_mag_m_actual+=stats::total_mag_m_actual;
		
		// increment data counter
		stats::data_counter+=1.0;

	}
	
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
	
	stats::total_mean_mag_m_actual=0.0;
	stats::total_mean_mag_m_norm=0.0;
	for(int mat=0;mat<mp::num_materials;mat++) stats::sublattice_mean_magm_array[mat]=0.0;
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
	
	// Reset Chi Data
	stats::MeanChi[0]=0.0;
	stats::MeanChi[1]=0.0;
	stats::MeanChi[2]=0.0;

	stats::MeanChiSquared[0]=0.0;
	stats::MeanChiSquared[1]=0.0;
	stats::MeanChiSquared[2]=0.0;

	stats::MeanChiDataCounter=0.0;
}

double max_torque(){
  //================================================================================================
  //
 //                                 subroutine torque
  //
  //                      Calculates total torque on the system
  //
  //================================================================================================

  //real(kind=dp) :: H (1:3)
  //real(kind=dp) :: S (1:3)
  
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
		MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&torque[0],3,MPI_DOUBLE,MPI_SUM);
		MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&stats::sublattice_mean_torque_x_array[0],mp::num_materials,MPI_DOUBLE,MPI_SUM);
		MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&stats::sublattice_mean_torque_y_array[0],mp::num_materials,MPI_DOUBLE,MPI_SUM);
		MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&stats::sublattice_mean_torque_z_array[0],mp::num_materials,MPI_DOUBLE,MPI_SUM);
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

///
/// @brief      Calculates the chi_l for the system
///
/// @return     void
///
void SystemSusceptibility(){

	double chi[3]={0.0,0.0,0.0};
	double chiM[3]={0.0,0.0,0.0};
	int chi_atoms=stats::num_atoms;
	
	for(int atom=0;atom<stats::num_atoms;atom++){

		// get atomic moment
		const int imat=atoms::type_array[atom];
		const double mu = mp::material[imat].mu_s_SI;
		
		chi[0]+=atoms::x_spin_array[atom];
		chi[1]+=atoms::y_spin_array[atom];
		chi[2]+=atoms::z_spin_array[atom];
		
	}

	// reduce torque on all nodes
	#ifdef MPICF
		MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&chi[0],3,MPI_DOUBLE,MPI_SUM);
		//MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&chisq[0],3,MPI_DOUBLE,MPI_SUM);
		MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&chi_atoms,1,MPI_INT,MPI_SUM);
	#endif

	stats::ChiAtoms=double(chi_atoms);
	double norm = 1.0/stats::ChiAtoms;
	
	chiM[0]=chi[0]*norm;
	chiM[1]=chi[1]*norm;
	chiM[2]=chi[2]*norm;
	
	// Calculate running totals
	stats::MeanChi[0]+=chiM[0];
	stats::MeanChi[1]+=chiM[1];
	stats::MeanChi[2]+=chiM[2];

	stats::MeanChiSquared[0]+=chiM[0]*chiM[0];
	stats::MeanChiSquared[1]+=chiM[1]*chiM[1];
	stats::MeanChiSquared[2]+=chiM[2]*chiM[2];

	stats::MeanChiDataCounter+=1.0;

	return;

}
	
} // End of Namespace
