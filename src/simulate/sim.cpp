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

   bool applied_field_enabled = true; // include the applied field in the spin field
   bool thermal_field_enabled = true; // include the thermal field in the spin field

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

   // the lists of atoms in each cell are only needed to initialise the dipole solvers; free them, as they hold an
   // inner vector for every cell in the system on every process
   std::vector < std::vector <int> >().swap(cells::index_atoms_array);
   std::vector < std::vector <double> >().swap(cells::atom_in_cell_coords_array_x);
   std::vector < std::vector <double> >().swap(cells::atom_in_cell_coords_array_y);
   std::vector < std::vector <double> >().swap(cells::atom_in_cell_coords_array_z);

   if(environment::enabled) environment::initialize(cs::system_dimensions[0],cs::system_dimensions[1],cs::system_dimensions[2]);

   // Precalculate initial statistics and then reset averages if not continuing a previous simulation. Done after the
   // dipole initialisation so that the initial magnetostatic energy includes the dipole field.
   // RE technically this double counts the last data point in the statistics, need to implement a reset_counter to fix.
   stats::update();
   if(!load_checkpoint_continue_flag) stats::reset();

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
   program::run();

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
