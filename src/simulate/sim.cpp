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
///
/// @file
/// @brief Contains the sim namespace and wrapper functions for system integration
///
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    09/03/2011
/// @internal
///   Created:    09/03/2011
///   Revision:     ---
///=====================================================================================
///

// Standard Libraries
#include <iostream>

// Vampire Header files
#include "anisotropy.hpp"
#include "atoms.hpp"
#include "program.hpp"
#include "cells.hpp"
#include "../cells/internal.hpp"
#include "../micromagnetic/internal.hpp"
#include "dipole.hpp"
#include "errors.hpp"
#include "gpu.hpp"
#include "grains.hpp"
#include "environment.hpp"
#include "hamr.hpp"
#include "integrate.hpp"
#include "material.hpp"
#include "montecarlo.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "spintorque.hpp"
#include "stats.hpp"
#include "stopwatch.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "vutil.hpp"
#include "micromagnetic.hpp"
#include "sld.hpp"

// sim module headers
#include "internal.hpp"


namespace sim{
   std::ofstream mag_file;

   int runs=1; /// for certain repetitions in programs

    //Global definition of some parameters in order to store them in chekcpoint files
   int64_t parity=-1;
   uint64_t output_atoms_file_counter=0;
   uint64_t output_cells_file_counter=0;
   uint64_t output_rate_counter=0;

   bool ext_demag=false;

   double Tmax=300.0;
   double Tmin=0.0;
   double Teq=300.0;
   double temperature=300.0;
   double delta_temperature=10.0;
   double H_applied=0.0;
   double H_vec[3]={0.0,0.0,1.0};
   double Hmin=-1.0; // T
   double Hmax=+1.0; // T
   double Hinc= 0.1; // T
   double Heq=0.0;
   double applied_field_angle_phi=0.0;
   double applied_field_angle_theta=0.0;
   bool applied_field_set_by_angle=false;

   double fmr_field_strength = 0.0; // Oscillating field strength (Tesla)
   double fmr_field_frequency = 1.0; // Oscillating field frequency (GHz)
   std::vector<double> fmr_field_unit_vector; // Oscillating field direction
   double fmr_field = 0.0; // Instantaneous value of the oscillating field strength H sin(wt)
   bool enable_fmr = false; // Flag to enable fmr field calculation

   double H=Hmax; // T
   int64_t iH=1; // uT
   // uint64_t iH=-1*vmath::iround(double(Hmax)*1.0E6); // uT

   double demag_factor[3]={0.0,0.0,0.0};
   // double head_position[2]={0.0,cs::system_dimensions[1]*0.5}; // A
   // double head_speed=30.0; /// nm/ns
   // bool   head_laser_on=false;
   bool   constraint_rotation=false; /// enables rotation of spins to new constraint direction
   bool   constraint_phi_changed=false; /// flag to note change in phi
   double constraint_phi=0.0; /// Constrained minimisation vector (azimuthal) [degrees]
   double constraint_phi_min=0.0; /// loop angle min [degrees]
   double constraint_phi_max=0.0; /// loop angle max [degrees]
   double constraint_phi_delta=5.0; /// loop angle delta [degrees]

   bool   constraint_theta_changed=false;
   double constraint_theta=0.0; /// Constrained minimisation vector (rotational) [degrees]
   double constraint_theta_min=0.0; /// loop angle min [degrees]
   double constraint_theta_max=0.0; // loop angle max [degrees]
   double constraint_theta_delta=5.0; /// loop angle delta [degrees]

   // LaGrange multiplier variables
   double lagrange_lambda_x=0.0;
   double lagrange_lambda_y=0.0;
   double lagrange_lambda_z=0.0;
   double lagrange_m=1.0;
   double lagrange_N=1000.0;
   bool   lagrange_multiplier=false;

   double cooling_time=100.0e-12; ///seconds
   int cooling_function_flag=0; /// 0 = exp, 1 = gaussian
   pump_functions_t pump_function=two_temperature;
   double pump_power=20.0; // mJ/cm^2;
   double pump_time=50.0e-15;
   double double_pump_power=20.0; // mJ/cm^2;
   double double_pump_Tmax=500.0;
   double double_pump_time=50.0e-15;
   double double_pump_delay=10.0e-12;
   double HeatSinkCouplingConstant=0.0; ///1.1e12 ~ sensible value
   double TTCe = 222.0; ///electron specific heat (gamma)
   double TTCl = 2.3E06; ///phonon specific heat
   double TTG = 6.6E17 ;///electron coupling constant
   double TTTe = 0.0; /// electron temperature
   double TTTp = 0.0; /// phonon temperature

   int system_simulation_flags;
   int hamiltonian_simulation_flags[10];

   bool local_temperature=false; /// flag to enable material specific temperature
   bool local_applied_field=false; /// flag to enable material specific applied field
   bool local_fmr_field=false; /// flag to enable material specific fmr field

   // Checkpoint flags and variables
   bool checkpoint_loaded_flag=false;  // Flag to determine if it is first step after loading checkpoint (true).
   bool load_checkpoint_flag=false; // Load spin configurations
   bool load_checkpoint_continue_flag=true; // Continue simulation from checkpoint time
   bool save_checkpoint_flag=false; // Save checkpoint
   bool save_checkpoint_continuous_flag=false; // save checkpoints during simulations
   int save_checkpoint_rate=1; // Default increment between checkpoints

   // Monte Carlo statistics counters
   double mc_statistics_moves = 0.0;
   double mc_statistics_reject = 0.0;

/// @brief Function to run one a single program
///
/// @callgraph
/// @callergraph
///
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.1
/// @date    09/03/2011
///
/// @return EXIT_SUCCESS
///
/// @internal
///   Created:    02/10/2008
///   Revision:   1.1 09/03/2011
///=====================================================================================
///
int run(){
   // Check for calling of function
   if(err::check==true) std::cout << "sim::run has been called" << std::endl;

   // Initialise simulation data structures
   sim::initialize(mp::num_materials);

   // Initialise integrate module data structures
   integrate::initialize();

   // Initialize vampire modules
   sim::internal::initialize_modules();

   montecarlo::initialize(atoms::num_atoms, grains::num_grains, atoms::grain_array);

   anisotropy::initialize(atoms::num_atoms, atoms::type_array, mp::mu_s_array);

   #ifdef MPICF
   if(sld::suzuki_trotter_parallel_initialized == false && sld::enabled) {
      sld::suzuki_trotter_parallel_init(atoms::x_coord_array, atoms::y_coord_array, atoms::z_coord_array,
                                        vmpi::min_dimensions, vmpi::max_dimensions);
   }
   #endif

   // now seed generator
   mtrandom::grnd.seed(vmpi::parallel_rng_seed(mtrandom::integration_seed));

   {
      // Set up statistical data sets
      #ifdef MPICF
         int num_atoms_for_statistics = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
      #else
         int num_atoms_for_statistics = atoms::num_atoms;
      #endif
      // unroll list of non-magnetic materials
      std::vector<bool> non_magnetic_materials_array(mp::num_materials, false);
      for(int m = 0; m < mp::num_materials; m++){
         if( mp::material[m].non_magnetic == 2 ) non_magnetic_materials_array[m] = true;
      }
      stats::initialize(num_atoms_for_statistics, mp::num_materials, grains::num_grains, atoms::m_spin_array, atoms::type_array, atoms::grain_array, atoms::category_array, non_magnetic_materials_array);
   }

   // Check for load spin configurations from checkpoint
   if(sim::load_checkpoint_flag) load_checkpoint();

   // Precalculate initial statistics and then reset averages if not continuing a previous simulation
   // RE technically this double counts the last data point in the statistics, need to implement a reset_counter to fix.
   stats::update();
   if(!load_checkpoint_continue_flag) stats::reset();

   // For continuous checkpoints inform user about I/O
   if(sim::save_checkpoint_continuous_flag) zlog << zTs() << "Continuously writing checkpoints to disk throughout simulation." << std::endl;

   // Initialize GPU acceleration if enabled
   if(gpu::acceleration) gpu::initialize();

   if (micromagnetic::discretisation_type > 0 || micromagnetic::internal::bias_magnets == true)
   micromagnetic::initialize(cells::num_local_cells,
                             cells::num_cells,
                             atoms::num_atoms,
                             mp::num_materials,
                             cells::atom_cell_id_array,
                             atoms::neighbour_list_array,
                             atoms::neighbour_list_start_index,
                             atoms::neighbour_list_end_index,
                             atoms::type_array,
                             mp::material,
                             atoms::x_coord_array,
                             atoms::y_coord_array,
                             atoms::z_coord_array,
                             cells::volume_array,
                             sim::temperature,
                             cells::num_atoms_in_unit_cell,
                             cs::system_dimensions[0],
                             cs::system_dimensions[1],
                             cs::system_dimensions[2],
                             cells::local_cell_array);

   // initialise dipole field calculation
   dipole::initialize(cells::num_atoms_in_unit_cell,
                     cells::num_cells,
                     cells::num_local_cells,
                     cells::macro_cell_size,
                     cells::local_cell_array,
                     cells::num_atoms_in_cell,
                     cells::num_atoms_in_cell_global, // <----
                     cells::index_atoms_array,
                     cells::volume_array,
                     cells::pos_and_mom_array,
                     cells::atom_in_cell_coords_array_x,
                     cells::atom_in_cell_coords_array_y,
                     cells::atom_in_cell_coords_array_z,
                     atoms::type_array,
                     cells::atom_cell_id_array,
                     atoms::x_coord_array,
                     atoms::y_coord_array,
                     atoms::z_coord_array,
                     atoms::x_spin_array,
                     atoms::y_spin_array,
                     atoms::z_spin_array,
                     atoms::m_spin_array,
                     atoms::magnetic,
                     atoms::num_atoms
   );

   if(environment::enabled) environment::initialize(cs::system_dimensions[0],cs::system_dimensions[1],cs::system_dimensions[2]);

   // For MPI version, calculate initialisation time
   if(vmpi::my_rank==0){
      #ifdef MPICF
         std::cout << "Time for initialisation: " << MPI_Wtime()-vmpi::start_time << std::endl;
         zlog << zTs() << "Time for initialisation: " << MPI_Wtime()-vmpi::start_time << std::endl;
         vmpi::start_time=MPI_Wtime(); // reset timer
      #endif
   }

   // Set timer for runtime
   stopwatch_t stopwatch;
   stopwatch.start();

   // Precondition spins at equilibration temperature
   montecarlo::monte_carlo_preconditioning();

   // For MPI version, calculate initialisation time
   if(vmpi::my_rank==0){
      std::cout << "Starting Simulation with Program ";
      zlog << zTs() << "Starting Simulation with Program ";
   }

   // Now set initial compute time
   #ifdef MPICF
   vmpi::ComputeTime=MPI_Wtime();
   vmpi::WaitTime=MPI_Wtime();
   vmpi::TotalComputeTime=0.0;
   vmpi::TotalWaitTime=0.0;
   #endif

   // Select program to run
   switch(program::program){
      case program::program_t::benchmark:
         if(vmpi::my_rank==0){
            std::cout << "Benchmark..." << std::endl;
            zlog << "Benchmark..." << std::endl;
         }
         program::bmark();
         break;

      case program::program_t::time_series:
         if(vmpi::my_rank==0){
            std::cout << "Time-Series..." << std::endl;
            zlog << "Time-Series..." << std::endl;
         }
         program::time_series();
         break;

      case program::program_t::hysteresis_loop:
         if(vmpi::my_rank==0){
            std::cout << "Hysteresis-Loop..." << std::endl;
            zlog << "Hysteresis-Loop..." << std::endl;
         }
         program::hysteresis();
         break;

      case program::program_t::static_hysteresis_loop:
         if(vmpi::my_rank==0){
            std::cout << "Static-Hysteresis-Loop..." << std::endl;
            zlog << "Static-Hysteresis-Loop..." << std::endl;
         }
         program::static_hysteresis();
         break;

      case program::program_t::curie_temperature:
         if(vmpi::my_rank==0){
            std::cout << "Curie-Temperature..." << std::endl;
            zlog << "Curie-Temperature..." << std::endl;
         }
         program::curie_temperature();
         break;

      case program::program_t::field_cool:
         if(vmpi::my_rank==0){
            std::cout << "Field-Cool..." << std::endl;
            zlog << "Field-Cool..." << std::endl;
         }
         program::field_cool();
         break;

      case program::program_t::laser_pulse:
         if(vmpi::my_rank==0){
            std::cout << "Temperature-Pulse..." << std::endl;
            zlog << "Temperature-Pulse..." << std::endl;
         }
         program::temperature_pulse();
         break;

      case program::program_t::hamr_simulation:
         if(vmpi::my_rank==0){
            std::cout << "HAMR-Simulation..." << std::endl;
            zlog << "HAMR-Simulation..." << std::endl;
         }
         program::hamr();
         break;

      case program::program_t::cmc_anisotropy:
         if(vmpi::my_rank==0){
            std::cout << "CMC-Anisotropy..." << std::endl;
            zlog << "CMC-Anisotropy..." << std::endl;
         }
         program::cmc_anisotropy();
         break;

      case program::program_t::hybrid_cmc:
         if(vmpi::my_rank==0){
            std::cout << "Hybrid-CMC..." << std::endl;
            zlog << "Hybrid-CMC..." << std::endl;
         }
         program::hybrid_cmc();
         break;

      case program::program_t::reverse_hybrid_cmc:
         if(vmpi::my_rank==0){
            std::cout << "Reverse-Hybrid-CMC..." << std::endl;
            zlog << "Reverse-Hybrid-CMC..." << std::endl;
         }
         program::reverse_hybrid_cmc();
         break;

      case program::program_t::lagrange_multiplier:
         if(vmpi::my_rank==0){
            std::cout << "LaGrange-Multiplier..." << std::endl;
            zlog << "LaGrange-Multiplier..." << std::endl;
         }
         program::lagrange_multiplier();
         break;

      case program::program_t::partial_hysteresis_loop:
         if(vmpi::my_rank==0){
            std::cout << "partial-hysteresis-loop..." << std::endl;
            zlog << "partial-hysteresis-loop..." << std::endl;
         }
         program::partial_hysteresis_loop();
         break;

      case program::program_t::localised_temperature_pulse:
         if(vmpi::my_rank==0){
            std::cout << "localised-temperature-pulse..." << std::endl;
            zlog << "localised-temperature-pulse..." << std::endl;
         }
         program::localised_temperature_pulse();
         break;

      case program::program_t::effective_damping:
         if(vmpi::my_rank==0){
            std::cout << "effective-damping..." << std::endl;
            zlog << "effective-damping..." << std::endl;
         }
         program::effective_damping();
         break;

      case program::program_t::fmr:
         if(vmpi::my_rank==0){
            std::cout << "fmr..." << std::endl;
            zlog << "fmr..." << std::endl;
         }
         program::fmr();
         break;

      case program::program_t::localised_field_cool:
         if(vmpi::my_rank==0){
            std::cout << "localised-field-cool..." << std::endl;
            zlog << "localised-field-cool..." << std::endl;
         }
         program::local_field_cool();
         break;

      case program::program_t::electrical_pulse:
         if(vmpi::my_rank==0){
            std::cout << "electrical-pulse..." << std::endl;
            zlog << "electrical-pulse..." << std::endl;
         }
         program::electrical_pulse();
         break;
      case program::program_t::field_pulse:
         if(vmpi::my_rank==0){
            std::cout << "field-pulse..." << std::endl;
            zlog << "field-pulse..." << std::endl;
         }
         program::field_pulse();
         break;
      case program::program_t::laser_electrical_pulse:
         if(vmpi::my_rank==0){
            std::cout << "laser-electrical-pulse..." << std::endl;
            zlog << "laser-electrical-pulse..." << std::endl;
         }
         program::laser_electrical_pulse();
         break;

      case program::program_t::diagnostic_boltzmann:
         if(vmpi::my_rank==0){
            std::cout << "Diagnostic-Boltzmann..." << std::endl;
            zlog << "Diagnostic-Boltzmann..." << std::endl;
         }
         program::boltzmann_dist();
         break;
      case program::program_t::setting:
         if(vmpi::my_rank==0){
            std::cout << "Setting..." << std::endl;
            zlog << "Setting..." << std::endl;
         }
         program::setting_process();
         break;
      //------------------------------------------------------------------------
      case program::program_t::domain_wall:
         if(vmpi::my_rank==0){
            std::cout << "Domain walls..." << std::endl;
            zlog << "Domain walls..." << std::endl;
         }
         program::domain_wall();
         break;
      //------------------------------------------------------------------------
      case program::program_t::exchange_stiffness:
         if(vmpi::my_rank==0){
            std::cout << "exchange stiffness..." << std::endl;
            zlog << "exchange stiffness..." << std::endl;
         }
         program::exchange_stiffness();
         break;
      //------------------------------------------------------------------------
      case program::program_t::mm_a_calculation:
         if(vmpi::my_rank==0){
            std::cout << "mm-A-calculation..." << std::endl;
            zlog << "mm-A-calculation..." << std::endl;
         }
         program::mm_A_calculation();
         break;
      //------------------------------------------------------------------------
      case program::program_t::field_sweep:
         if(vmpi::my_rank==0){
            std::cout << "field-sweep..." << std::endl;
            zlog << "field-sweep..." << std::endl;
         }
         program::field_sweep();
      break;
      //------------------------------------------------------------------------
      case program::program_t::disk_tracks:
         if(vmpi::my_rank==0){
            std::cout << "Tracks..." << std::endl;
            zlog << "Tracks..." << std::endl;
         }
         program::tracks();
         break;
      //------------------------------------------------------------------------
      case program::program_t::diagnostic_boltzmann_micromagnetic_llg:
         if(vmpi::my_rank==0){
            std::cout << "diagnostic-boltzmann-micromganetic-llg..." << std::endl;
            zlog << "diagnostic-boltzmann-micromganetic-llg..." << std::endl;
         }
         program::boltzmann_dist_micromagnetic_llg();
         break;
      // JRH call spin_waves program -------------------------------------------
      case program::program_t::spin_waves:
         if(vmpi::my_rank==0){
            std::cout << "Spin-waves..." << std::endl;
            zlog << "Spin-waves..." << std::endl;
         }
         program::spin_waves();
         break;
      //------------------------------------------------------------------------
      case program::program_t::hamr_cool:
         if(vmpi::my_rank==0){
            std::cout << "HAMR cool..." << std::endl;
            zlog << "HAMR cool..." << std::endl;
         }
         program::hamr_cool();
         break;

      //------------------------------------------------------------------------
      default:{
         std::cerr << "Unknown Internal Program ID "<< static_cast<int>(program::program) << " requested, exiting" << std::endl;
         zlog << "Unknown Internal Program ID "<< static_cast<int>(program::program) << " requested, exiting" << std::endl;
         err::vexit();
         }
   }

   std::cout <<     "Simulation run time [s]: " << stopwatch.elapsed_seconds() << std::endl;
   zlog << zTs() << "Simulation run time [s]: " << stopwatch.elapsed_seconds() << std::endl;

   //------------------------------------------------
   // Output Monte Carlo statistics if applicable
   //------------------------------------------------
   if(integrate::integrator == integrate::integrator_t::monte_carlo){
      std::cout << "Monte Carlo statistics:" << std::endl;
      std::cout << "\tTotal moves: " << long(sim::mc_statistics_moves) << std::endl;
      std::cout << "\t" << ((sim::mc_statistics_moves - sim::mc_statistics_reject)/sim::mc_statistics_moves)*100.0 << "% Accepted" << std::endl;
      std::cout << "\t" << (sim::mc_statistics_reject/sim::mc_statistics_moves)*100.0                              << "% Rejected" << std::endl;
      zlog << zTs() << "Monte Carlo statistics:" << std::endl;
      zlog << zTs() << "\tTotal moves: " << sim::mc_statistics_moves << std::endl;
      zlog << zTs() << "\t" << ((sim::mc_statistics_moves - sim::mc_statistics_reject)/sim::mc_statistics_moves)*100.0 << "% Accepted" << std::endl;
      zlog << zTs() << "\t" << (sim::mc_statistics_reject/sim::mc_statistics_moves)*100.0                              << "% Rejected" << std::endl;
   }
   if(integrate::integrator == integrate::integrator_t::cmc || integrate::integrator == integrate::integrator_t::hybrid_cmc){
      std::cout << "Constrained Monte Carlo statistics:" << std::endl;
      std::cout << "\tTotal moves: " << montecarlo::cmc::mc_total << std::endl;
      std::cout << "\t" << (montecarlo::cmc::mc_success/montecarlo::cmc::mc_total)*100.0    << "% Accepted" << std::endl;
      std::cout << "\t" << (montecarlo::cmc::energy_reject/montecarlo::cmc::mc_total)*100.0 << "% Rejected (Energy)" << std::endl;
      std::cout << "\t" << (montecarlo::cmc::sphere_reject/montecarlo::cmc::mc_total)*100.0 << "% Rejected (Sphere)" << std::endl;
      zlog << zTs() << "Constrained Monte Carlo statistics:" << std::endl;
      zlog << zTs() << "\tTotal moves: " << montecarlo::cmc::mc_total << std::endl;
      zlog << zTs() << "\t" << (montecarlo::cmc::mc_success/montecarlo::cmc::mc_total)*100.0    << "% Accepted" << std::endl;
      zlog << zTs() << "\t" << (montecarlo::cmc::energy_reject/montecarlo::cmc::mc_total)*100.0 << "% Rejected (Energy)" << std::endl;
      zlog << zTs() << "\t" << (montecarlo::cmc::sphere_reject/montecarlo::cmc::mc_total)*100.0 << "% Rejected (Sphere)" << std::endl;
   }

   //program::LLB_Boltzmann();

   // De-initialize GPU
   if(gpu::acceleration) gpu::finalize();

   // optionally save checkpoint file
   if(sim::save_checkpoint_flag && !sim::save_checkpoint_continuous_flag) save_checkpoint();

   return EXIT_SUCCESS;
}

/// @brief Wrapper function to call integrators
///
/// @callgraph
/// @callergraph
///
/// @details Calls serial or parallel integrator
///
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    05/02/2011
///
/// @return EXIT_SUCCESS
///
/// @internal
///   Created:    05/02/2011
///   Revision:     ---
///=====================================================================================
///
int integrate(uint64_t n_steps){

   // Check for calling of function
   if(err::check==true) std::cout << "sim::integrate has been called" << std::endl;

   // Call integrate module
   return integrate::integrate(n_steps);
}

} // Namespace sim
