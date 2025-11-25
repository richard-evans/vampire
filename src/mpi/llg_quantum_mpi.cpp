//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) David R Papp 2024. All rights reserved.
//
//------------------------------------------------------------------------------
//

#ifdef MPICF
// Standard Libraries
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <complex>
#include <iomanip>
#include <numeric>

// Library for FFT
#ifdef FFT
#include <fftw3.h>
#endif

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "LLG.hpp"
#include "material.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vmpi.hpp"
#include "quantum.hpp"

#include "../simulate/internal.hpp"




int calculate_spin_fields(const int,const int);
int calculate_external_fields(const int,const int);

namespace sim{


    int LLGQ_mpi_init(){
      // check calling of routine if error checking is activated
      if(err::check==true){std::cout << "sim::LLG_init has been called" << std::endl;}
      using namespace LLGQ_arrays;


      x_w_array.resize(atoms::num_atoms, 0.0);
      y_w_array.resize(atoms::num_atoms, 0.0);
      z_w_array.resize(atoms::num_atoms, 0.0);
      x_v_array.resize(atoms::num_atoms, 0.0);
      y_v_array.resize(atoms::num_atoms, 0.0);
      z_v_array.resize(atoms::num_atoms, 0.0);

      k1_storage.resize(atoms::num_atoms, std::vector<double>(9));
      k2_storage.resize(atoms::num_atoms, std::vector<double>(9));
      k3_storage.resize(atoms::num_atoms, std::vector<double>(9));
      k4_storage.resize(atoms::num_atoms, std::vector<double>(9));
      y_pred_storage.resize(atoms::num_atoms, std::vector<double>(9));
      y_in_storage.resize(atoms::num_atoms, std::vector<double>(9));


      // Initialize Quantum Module
      quantum::initialize();

      LLG_set=true;

      
      vmpi::barrier();

      return EXIT_SUCCESS;
   }

   namespace internal{
    void spinDynamics(const double* y, const double* H, double* dydt);
   }


    int llg_quantum_mpi_step(){
        //======================================================
        // Subroutine to perform a single LSF integration step
        //======================================================


        //----------------------------------------------------------
        // check calling of routine if error checking is activated
        //----------------------------------------------------------
        if(err::check==true){std::cout << "LLG_quantum_mpi has been called" << std::endl;}

        using namespace LLGQ_arrays;
        using namespace sim;
        using namespace sim::internal;
        // Check for initialisation of LSF integration arrays
        if(LLG_set==false) sim::LLGQ_mpi_init();
	
	
        
        //----------------------------------------
        // Local variables for system generation
        //----------------------------------------
        //const int num_atoms = atoms::num_atoms;
        const int pre_comm_si = 0;
        const int pre_comm_ei = vmpi::num_core_atoms;
        const int post_comm_si = vmpi::num_core_atoms;
        const int post_comm_ei = vmpi::num_core_atoms+vmpi::num_bdry_atoms;

        const double dt = mp::dt;
        const double half_dt = 0.5 * mp::dt;
        const double dt_over_6 = mp::dt / 6.0;
        const int M = M_decimation;

        std::vector<double> H(3);

        
        //----------------------------------------
		// Initiate halo swap
		//----------------------------------------
		
        vmpi::mpi_init_halo_swap();

        //----------------------------------------
		// Calculate fields (core)
		//----------------------------------------

		calculate_spin_fields(pre_comm_si,pre_comm_ei);
        calculate_external_fields(pre_comm_si,post_comm_ei);

        //----------------------------------------
        // Store initial spin positions (all)
        //----------------------------------------

        for(int atom=pre_comm_si;atom<post_comm_ei;atom++){
            y_in_storage[atom][0] = atoms::x_spin_array[atom];
            y_in_storage[atom][1] = atoms::y_spin_array[atom];
            y_in_storage[atom][2] = atoms::z_spin_array[atom];
            y_in_storage[atom][3] = LLGQ_arrays::x_v_array[atom];
            y_in_storage[atom][4] = LLGQ_arrays::y_v_array[atom];
            y_in_storage[atom][5] = LLGQ_arrays::z_v_array[atom];
            y_in_storage[atom][6] = LLGQ_arrays::x_w_array[atom];
            y_in_storage[atom][7] = LLGQ_arrays::y_w_array[atom];
            y_in_storage[atom][8] = LLGQ_arrays::z_w_array[atom];
        }

        //----------------------------------------
        // Calculate K1 (Core)
        //----------------------------------------

        for(int atom=pre_comm_si;atom<pre_comm_ei;atom++){
            H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + quantum::get_field(atom, 0, 0.0);
            H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + quantum::get_field(atom, 1, 0.0);
            H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + quantum::get_field(atom, 2, 0.0);

            // Calculate K1
            spinDynamics(y_in_storage[atom].data(), H.data(), k1_storage[atom].data());

            // Calculate y_pred for k2
            for (size_t i = 0; i < 9; ++i) {
               y_pred_storage[atom][i] = y_in_storage[atom][i] + half_dt * k1_storage[atom][i];
            }

            // Normalize spin length
            double S_magnitude = std::sqrt(y_pred_storage[atom][0] * y_pred_storage[atom][0] + y_pred_storage[atom][1] * y_pred_storage[atom][1] + y_pred_storage[atom][2] * y_pred_storage[atom][2]);
            const double inv_mag = 1.0 / S_magnitude;
            y_pred_storage[atom][0] *= inv_mag;
            y_pred_storage[atom][1] *= inv_mag;
            y_pred_storage[atom][2] *= inv_mag;
        
        }

        //----------------------------------------
        // Complete halo swap
        //----------------------------------------
        
        vmpi::mpi_complete_halo_swap();

        //----------------------------------------
		// Calculate fields (boundary)
		//----------------------------------------

		calculate_spin_fields(post_comm_si,post_comm_ei);
		calculate_external_fields(post_comm_si,post_comm_ei);

        //----------------------------------------
        // Calculate K1 (boundary)
        //----------------------------------------

        for(int atom=post_comm_si;atom<post_comm_ei;atom++){
            H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + quantum::get_field(atom, 0, 0.0);
            H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + quantum::get_field(atom, 1, 0.0);
            H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + quantum::get_field(atom, 2, 0.0);

            // Calculate K1
            spinDynamics(y_in_storage[atom].data(), H.data(), k1_storage[atom].data());

            // Calculate y_pred for k2
            for (size_t i = 0; i < 9; ++i) {
               y_pred_storage[atom][i] = y_in_storage[atom][i] + half_dt * k1_storage[atom][i];
            }

            // Normalize spin length
            double S_magnitude = std::sqrt(y_pred_storage[atom][0] * y_pred_storage[atom][0] + y_pred_storage[atom][1] * y_pred_storage[atom][1] + y_pred_storage[atom][2] * y_pred_storage[atom][2]);
            const double inv_mag = 1.0 / S_magnitude;
            y_pred_storage[atom][0] *= inv_mag;
            y_pred_storage[atom][1] *= inv_mag;
            y_pred_storage[atom][2] *= inv_mag;
        }

        //------------------------------------------
        // Update spin positions (all)
        //------------------------------------------

        for(int atom=pre_comm_si;atom<post_comm_ei;atom++){
            // Update spin for field calculation
            atoms::x_spin_array[atom] = y_pred_storage[atom][0];
            atoms::y_spin_array[atom] = y_pred_storage[atom][1];
            atoms::z_spin_array[atom] = y_pred_storage[atom][2];
        }

        //------------------------------------------
        // Initiate second halo swap
        //------------------------------------------
        vmpi::mpi_init_halo_swap();

        //------------------------------------------
        // Recalculate spin dependent fields (core)
        //------------------------------------------

        calculate_spin_fields(pre_comm_si,pre_comm_ei);
        
        //----------------------------------------
        // Calculate K2 (core)
        //----------------------------------------

        for(int atom=pre_comm_si;atom<pre_comm_ei;atom++){
            H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + quantum::get_field(atom, 0, 0.5);
            H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + quantum::get_field(atom, 1, 0.5);
            H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + quantum::get_field(atom, 2, 0.5);

            // Calculate K2
            spinDynamics(y_pred_storage[atom].data(), H.data(), k2_storage[atom].data());

            // Calculate y_pred for k3
            for (size_t i = 0; i < 9; ++i) {
               y_pred_storage[atom][i] = y_in_storage[atom][i] + half_dt * k2_storage[atom][i];
            }

            // Normalize spin length
            double S_magnitude = std::sqrt(y_pred_storage[atom][0] * y_pred_storage[atom][0] + y_pred_storage[atom][1] * y_pred_storage[atom][1] + y_pred_storage[atom][2] * y_pred_storage[atom][2]);
            const double inv_mag = 1.0 / S_magnitude;
            y_pred_storage[atom][0] *= inv_mag;
            y_pred_storage[atom][1] *= inv_mag;
            y_pred_storage[atom][2] *= inv_mag;
        }

        //------------------------------------------
        // Complete second halo swap
        //------------------------------------------
        
        vmpi::mpi_complete_halo_swap();

        //------------------------------------------
		// Recalculate spin dependent fields (boundary)
		//------------------------------------------

		calculate_spin_fields(post_comm_si,post_comm_ei);
    
        //----------------------------------------
        // Calculate K2 (boundary)
        //----------------------------------------

        for(int atom=post_comm_si;atom<post_comm_ei;atom++){
            H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + quantum::get_field(atom, 0, 0.5);
            H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + quantum::get_field(atom, 1, 0.5);
            H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + quantum::get_field(atom, 2, 0.5);

            // Calculate K2
            spinDynamics(y_pred_storage[atom].data(), H.data(), k2_storage[atom].data());

            // Calculate y_pred for k3
            for (size_t i = 0; i < 9; ++i) {
               y_pred_storage[atom][i] = y_in_storage[atom][i] + half_dt * k2_storage[atom][i];
            }

            // Normalize spin length
            double S_magnitude = std::sqrt(y_pred_storage[atom][0] * y_pred_storage[atom][0] + y_pred_storage[atom][1] * y_pred_storage[atom][1] + y_pred_storage[atom][2] * y_pred_storage[atom][2]);
            const double inv_mag = 1.0 / S_magnitude;
            y_pred_storage[atom][0] *= inv_mag;
            y_pred_storage[atom][1] *= inv_mag;
            y_pred_storage[atom][2] *= inv_mag;
        }

        //------------------------------------------
        // Update spin positions (all)
        //------------------------------------------

        for(int atom=pre_comm_si;atom<post_comm_ei;atom++){
            // Update spin for field calculation
            atoms::x_spin_array[atom] = y_pred_storage[atom][0];
            atoms::y_spin_array[atom] = y_pred_storage[atom][1];
            atoms::z_spin_array[atom] = y_pred_storage[atom][2];
        }

        //------------------------------------------
        // Initiate third halo swap
        //------------------------------------------
        vmpi::mpi_init_halo_swap();

        //------------------------------------------
        // Recalculate spin dependent fields (core)
        //------------------------------------------

        calculate_spin_fields(pre_comm_si,pre_comm_ei);

        //----------------------------------------
        // Calculate K3 (core)
        //----------------------------------------

        for(int atom=pre_comm_si;atom<pre_comm_ei;atom++){
            H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + quantum::get_field(atom, 0, 0.5);
            H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + quantum::get_field(atom, 1, 0.5);
            H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + quantum::get_field(atom, 2, 0.5);

            // Calculate K3
            spinDynamics(y_pred_storage[atom].data(), H.data(), k3_storage[atom].data());

            // Calculate y_pred for k4
            for (size_t i = 0; i < 9; ++i) {
               y_pred_storage[atom][i] = y_in_storage[atom][i] + dt * k3_storage[atom][i];
            }

            // Normalize spin length
            double S_magnitude = std::sqrt(y_pred_storage[atom][0] * y_pred_storage[atom][0] + y_pred_storage[atom][1] * y_pred_storage[atom][1] + y_pred_storage[atom][2] * y_pred_storage[atom][2]);
            const double inv_mag = 1.0 / S_magnitude;
            y_pred_storage[atom][0] *= inv_mag;
            y_pred_storage[atom][1] *= inv_mag;
            y_pred_storage[atom][2] *= inv_mag;
        }

        //------------------------------------------
        // Complete third halo swap
        //------------------------------------------
        
        vmpi::mpi_complete_halo_swap();
;
        //------------------------------------------
        // Recalculate spin dependent fields (boundary)
        //------------------------------------------
        
        calculate_spin_fields(post_comm_si,post_comm_ei);

        //----------------------------------------
        // Calculate K3 (boundary)
        //----------------------------------------

        for(int atom=post_comm_si;atom<post_comm_ei;atom++){
            H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + quantum::get_field(atom, 0, 0.5);
            H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + quantum::get_field(atom, 1, 0.5);
            H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + quantum::get_field(atom, 2, 0.5);

            // Calculate K3
            spinDynamics(y_pred_storage[atom].data(), H.data(), k3_storage[atom].data());

            // Calculate y_pred for k4
            for (size_t i = 0; i < 9; ++i) {
               y_pred_storage[atom][i] = y_in_storage[atom][i] + dt * k3_storage[atom][i];
            }

            // Normalize spin length
            double S_magnitude = std::sqrt(y_pred_storage[atom][0] * y_pred_storage[atom][0] + y_pred_storage[atom][1] * y_pred_storage[atom][1] + y_pred_storage[atom][2] * y_pred_storage[atom][2]);
            const double inv_mag = 1.0 / S_magnitude;
            y_pred_storage[atom][0] *= inv_mag;
            y_pred_storage[atom][1] *= inv_mag;
            y_pred_storage[atom][2] *= inv_mag;
        }

        //------------------------------------------
        // Update spin positions (all)
        //------------------------------------------

        for(int atom=pre_comm_si;atom<post_comm_ei;atom++){
            // Update spin for field calculation
            atoms::x_spin_array[atom] = y_pred_storage[atom][0];
            atoms::y_spin_array[atom] = y_pred_storage[atom][1];
            atoms::z_spin_array[atom] = y_pred_storage[atom][2];
        }

        //------------------------------------------
        // Initiate fourth halo swap
        //------------------------------------------
        vmpi::mpi_init_halo_swap();

        //------------------------------------------
		// Recalculate spin dependent fields (core)
		//------------------------------------------

		calculate_spin_fields(pre_comm_si,pre_comm_ei);

        //----------------------------------------
        // Calculate K4 (core)
        //----------------------------------------

        for(int atom=pre_comm_si;atom<pre_comm_ei;atom++){
            H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + quantum::get_field(atom, 0, 1.0);
            H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + quantum::get_field(atom, 1, 1.0);
            H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + quantum::get_field(atom, 2, 1.0);

            // Calculate K4
            spinDynamics(y_pred_storage[atom].data(), H.data(), k4_storage[atom].data());
        }

        //------------------------------------------
        // Complete fourth halo swap
        //------------------------------------------
        
        vmpi::mpi_complete_halo_swap();

        //------------------------------------------
        // Recalculate spin dependent fields (boundary)
        //------------------------------------------
		
        calculate_spin_fields(post_comm_si,post_comm_ei);

        //----------------------------------------
        // Calculate K4 (boundary)
        //----------------------------------------

        for(int atom=post_comm_si;atom<post_comm_ei;atom++){
            H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + quantum::get_field(atom, 0, 1.0);
            H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + quantum::get_field(atom, 1, 1.0);
            H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + quantum::get_field(atom, 2, 1.0);

            // Calculate K4
            spinDynamics(y_pred_storage[atom].data(), H.data(), k4_storage[atom].data());
        }

        //----------------------------------------
        // Calculate RK4 Step
        //----------------------------------------

        for(int atom=pre_comm_si;atom<post_comm_ei;atom++){
            // Final update and normalization for spin components
            for (size_t i = 0; i < 3; ++i) {
               y_pred_storage[atom][i] = y_in_storage[atom][i] + dt_over_6 * (k1_storage[atom][i] + 2.0 * k2_storage[atom][i] + 2.0 * k3_storage[atom][i] + k4_storage[atom][i]);
            }

            double S_magnitude = std::sqrt(y_pred_storage[atom][0] * y_pred_storage[atom][0] + y_pred_storage[atom][1] * y_pred_storage[atom][1] + y_pred_storage[atom][2] * y_pred_storage[atom][2]);
            const double inv_mag = 1.0 / S_magnitude;
            y_pred_storage[atom][0] *= inv_mag;
            y_pred_storage[atom][1] *= inv_mag;
            y_pred_storage[atom][2] *= inv_mag;

            atoms::x_spin_array[atom] = y_pred_storage[atom][0];
            atoms::y_spin_array[atom] = y_pred_storage[atom][1];
            atoms::z_spin_array[atom] = y_pred_storage[atom][2];

            // Update auxiliary variables
            LLGQ_arrays::x_v_array[atom] = y_in_storage[atom][3] + (dt_over_6) * (k1_storage[atom][3] + 2.0 * k2_storage[atom][3] + 2.0 * k3_storage[atom][3] + k4_storage[atom][3]);
            LLGQ_arrays::y_v_array[atom] = y_in_storage[atom][4] + (dt_over_6) * (k1_storage[atom][4] + 2.0 * k2_storage[atom][4] + 2.0 * k3_storage[atom][4] + k4_storage[atom][4]);
            LLGQ_arrays::z_v_array[atom] = y_in_storage[atom][5] + (dt_over_6) * (k1_storage[atom][5] + 2.0 * k2_storage[atom][5] + 2.0 * k3_storage[atom][5] + k4_storage[atom][5]);
            LLGQ_arrays::x_w_array[atom] = y_in_storage[atom][6] + (dt_over_6) * (k1_storage[atom][6] + 2.0 * k2_storage[atom][6] + 2.0 * k3_storage[atom][6] + k4_storage[atom][6]);
            LLGQ_arrays::y_w_array[atom] = y_in_storage[atom][7] + (dt_over_6) * (k1_storage[atom][7] + 2.0 * k2_storage[atom][7] + 2.0 * k3_storage[atom][7] + k4_storage[atom][7]);
            LLGQ_arrays::z_w_array[atom] = y_in_storage[atom][8] + (dt_over_6) * (k1_storage[atom][8] + 2.0 * k2_storage[atom][8] + 2.0 * k3_storage[atom][8] + k4_storage[atom][8]);
        }

        // Swap timers compute -> wait
        vmpi::TotalComputeTime+=vmpi::SwapTimer(vmpi::ComputeTime, vmpi::WaitTime);

        // Wait for other processors
        vmpi::barrier();

        // Increment noise index
        quantum::increment_time();

        // Swap timers wait -> compute
        vmpi::TotalWaitTime+=vmpi::SwapTimer(vmpi::WaitTime, vmpi::ComputeTime);

    return EXIT_SUCCESS;
    }
} // end of namespace sim
#endif
