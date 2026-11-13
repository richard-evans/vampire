//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2014. All rights reserved.
//
//------------------------------------------------------------------------------

// System headers
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

// Program headers
#include "atoms.hpp"
#include "create.hpp"
#include "errors.hpp"
#include "gpu.hpp"
#include "hamr.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "sld.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "vutil.hpp"
#include "montecarlo.hpp"
#include "program.hpp"

namespace{

   // Place a magic key at the end of the file to indicate if the file has atom keys stored
   const uint64_t keyed_section_magic = 0x31594b5043504d56ULL;

   std::string checkpoint_file_name(const int rank){
      std::stringstream chkfilenamess;
      chkfilenamess << "vampire" << rank << ".chk";
      return chkfilenamess.str();
   }

   //--------------------------------------------------------------------------
   // Function to return the position of the section with the atom keys of a
   // checkpoint file, or 0 if the file has none (older files)
   //--------------------------------------------------------------------------
   uint64_t keyed_section_position(std::ifstream& chkfile){
      chkfile.clear();
      chkfile.seekg(0, std::ios::end);
      const std::streamoff size = chkfile.tellg();
      if(size < std::streamoff(2*sizeof(uint64_t))) return 0;
      uint64_t position = 0;
      uint64_t magic = 0;
      chkfile.seekg(size - std::streamoff(2*sizeof(uint64_t)));
      chkfile.read((char*)&position,sizeof(uint64_t));
      chkfile.read((char*)&magic,sizeof(uint64_t));
      return (chkfile && magic == keyed_section_magic) ? position : 0;
   }

   //--------------------------------------------------------------------------
   // Function to load the spins of the local atoms from the checkpoint files
   // of another decomposition: each atom is found by its key in the files of
   // the processes whose region overlaps the region of this process
   //--------------------------------------------------------------------------
   void load_keyed_spins(){

      const uint64_t num_local_atoms = uint64_t(atoms::num_atoms - vmpi::num_halo_atoms);

      // local atoms sorted by key
      std::vector<std::pair<uint64_t,uint64_t> > order(num_local_atoms);
      for(uint64_t atom = 0; atom < num_local_atoms; atom++) order[atom] = std::make_pair(create::atom_keys[atom], atom);
      std::sort(order.begin(), order.end());

      // only the regions of a geometric decomposition limit the files to read
      const bool own_region = vmpi::num_processors > 1 && vmpi::mpi_mode == 0;
      const double tolerance = 1.0e-3; // Angstroms

      std::vector<int> found(num_local_atoms, 0);
      uint64_t wrong_material = 0;
      uint64_t atoms_in_files = 0;

      // the checkpoint consists of the files of the processes that wrote file 0
      // (prevents confusion from other files that may be left over from earlier simulations)
      int64_t num_files = 1;
      for(int64_t file = 0; file < num_files; file++){

         const std::string chkfilename = checkpoint_file_name(file);
         std::ifstream chkfile(chkfilename.c_str(), std::ios::binary);
         if(!chkfile.is_open()) err::fatal("Error: Unable to open checkpoint file " + chkfilename + " for reading. Exiting.");
         const uint64_t section = keyed_section_position(chkfile);
         if(section == 0) err::fatal("Error: Checkpoint file " + chkfilename + " holds no atom keys and cannot be loaded into a different decomposition. Exiting.");

         int64_t num_processors64 = 0;
         double min_dimensions[3];
         double max_dimensions[3];
         uint64_t spins_position = 0;
         uint64_t natoms64 = 0;
         chkfile.seekg(section);
         chkfile.read((char*)&num_processors64,sizeof(int64_t));
         chkfile.read((char*)&min_dimensions[0],sizeof(double)*3);
         chkfile.read((char*)&max_dimensions[0],sizeof(double)*3);
         chkfile.read((char*)&spins_position,sizeof(uint64_t));
         chkfile.read((char*)&natoms64,sizeof(uint64_t));
         if(file == 0) num_files = num_processors64;
         if(num_processors64 != num_files){
            err::fatal("Error: Checkpoint file " + chkfilename + " was written by " + std::to_string(num_processors64) + " processes and " +
                       checkpoint_file_name(0) + " by " + std::to_string(num_files) + ": they belong to different simulations. Exiting.");
         }
         atoms_in_files += natoms64;

         if(own_region && num_processors64 > 1){
            bool overlap = true;
            for(int i = 0; i < 3; i++){
               overlap = overlap && min_dimensions[i] <= vmpi::max_dimensions[i] + tolerance && vmpi::min_dimensions[i] <= max_dimensions[i] + tolerance;
            }
            if(!overlap) continue;
         }

         std::vector<uint64_t> keys(natoms64);
         std::vector<int32_t> materials(natoms64);
         std::vector<double> sx(natoms64), sy(natoms64), sz(natoms64);
         chkfile.read((char*)keys.data(),sizeof(uint64_t)*natoms64);
         chkfile.read((char*)materials.data(),sizeof(int32_t)*natoms64);
         chkfile.seekg(spins_position);
         chkfile.read((char*)sx.data(),sizeof(double)*natoms64);
         chkfile.read((char*)sy.data(),sizeof(double)*natoms64);
         chkfile.read((char*)sz.data(),sizeof(double)*natoms64);
         if(!chkfile) err::fatal("Error: Checkpoint file " + chkfilename + " is incomplete. Exiting.");

         for(uint64_t i = 0; i < natoms64; i++){
            const auto it = std::lower_bound(order.begin(), order.end(), std::make_pair(keys[i], uint64_t(0)));
            if(it == order.end() || it->first != keys[i]) continue; // atom of another process
            const uint64_t atom = it->second;
            found[atom]++;
            if(materials[i] != atoms::type_array[atom]) wrong_material++;
            atoms::x_spin_array[atom] = sx[i];
            atoms::y_spin_array[atom] = sy[i];
            atoms::z_spin_array[atom] = sz[i];
         }

      }

      // every local atom must be found exactly once, with the same material, and the files
      // must hold no other atoms
      uint64_t errors[4] = {0, 0, wrong_material, num_local_atoms};
      for(uint64_t atom = 0; atom < num_local_atoms; atom++){
         if(found[atom] == 0) errors[0]++;
         if(found[atom] > 1) errors[1]++;
      }
      #ifdef MPICF
         MPI_Allreduce(MPI_IN_PLACE, errors, 4, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
      #endif
      if(errors[0] + errors[1] + errors[2] > 0 || atoms_in_files != errors[3]){
         std::ostringstream message;
         message << "Error: The atoms of the checkpoint files (" << atoms_in_files << ") differ from the generated atoms (" << errors[3] << "): "
                 << errors[0] << " atoms missing, " << errors[1] << " atoms found more than once, " << errors[2] << " atoms of a different material. Random choices of "
                 << "creation (alloying, intermixing, dilution) differ between numbers of processes unless "
                 << "create:decomposition-independent-random is set when writing and loading the checkpoint. Exiting.";
         err::fatal(message.str());
      }

   }

}

//-----------------------------------------------------------------------------
// Function to save checkpoint file
//-----------------------------------------------------------------------------
void save_checkpoint(){

   // For gpu acceleration get statistics from device
   if(gpu::acceleration) gpu::stats::get();

   // convert number of atoms, rank and time to standard long int
   uint64_t natoms64 = uint64_t(atoms::num_atoms-vmpi::num_halo_atoms);
   int64_t time64 = int64_t(sim::time);
   int64_t eqtime64 = int64_t(sim::equilibration_time);
   int64_t parity64 = int64_t(sim::parity);
   int64_t iH64 = int64_t(sim::iH);
   double temp = sim::temperature;
   int64_t output_atoms_file_counter64 = int64_t(sim::output_atoms_file_counter);
   int64_t output_cells_file_counter64 = int64_t(sim::output_cells_file_counter);
   int64_t output_rate_counter64 = int64_t(sim::output_rate_counter);
   double constr_theta = montecarlo::constraint_theta;
   double constr_phi   = montecarlo::constraint_phi;
   bool flag_constraint_theta_changed = montecarlo::constraint_theta_changed;
   bool flag_constraint_phi_changed   = montecarlo::constraint_phi_changed;

   // determine checkpoint file name
   std::stringstream chkfilenamess;
   chkfilenamess << "vampire" << vmpi::my_rank << ".chk";
   std::string chkfilename = chkfilenamess.str();

   // open checkpoint file
   std::ofstream chkfile;
   chkfile.open(chkfilename.c_str(),std::ios::binary);

   // check for open file
   if(!chkfile.is_open()){
      std::ostringstream message;
      message << "Error: Unable to open checkpoint file " << chkfilename << " for writing. Exiting.";
      err::fatal(message.str());
   }

   // get state of random number generator
   std::vector<uint32_t> mt_state(624); // 624 is hard coded in mt implementation. uint64 assumes same size as unsigned long
   int32_t mt_p=0; // position in rng state
   mt_p=mtrandom::grnd.get_state(mt_state);

   // write checkpoint variables to file
   chkfile.write(reinterpret_cast<const char*>(&natoms64),sizeof(uint64_t));
   chkfile.write(reinterpret_cast<const char*>(&time64),sizeof(int64_t));
   chkfile.write(reinterpret_cast<const char*>(&eqtime64),sizeof(int64_t));
   chkfile.write(reinterpret_cast<const char*>(&parity64),sizeof(int64_t));
   chkfile.write(reinterpret_cast<const char*>(&iH64),sizeof(int64_t));
   chkfile.write(reinterpret_cast<const char*>(&temp),sizeof(double));
   chkfile.write(reinterpret_cast<const char*>(&constr_theta),sizeof(double));
   chkfile.write(reinterpret_cast<const char*>(&constr_phi),sizeof(double));
   chkfile.write(reinterpret_cast<const char*>(&flag_constraint_theta_changed),sizeof(bool));
   chkfile.write(reinterpret_cast<const char*>(&flag_constraint_phi_changed  ),sizeof(bool));
   chkfile.write(reinterpret_cast<const char*>(&output_atoms_file_counter64),sizeof(int64_t));
   chkfile.write(reinterpret_cast<const char*>(&output_cells_file_counter64),sizeof(int64_t));
   chkfile.write(reinterpret_cast<const char*>(&output_rate_counter64),sizeof(int64_t));
   chkfile.write(reinterpret_cast<const char*>(&mt_p),sizeof(int32_t));
   chkfile.write(reinterpret_cast<const char*>(&mt_state[0]),sizeof(uint32_t)*mt_state.size());

   // write spin array to file
   const uint64_t spins_position = uint64_t(chkfile.tellp());
   chkfile.write(reinterpret_cast<const char*>(&atoms::x_spin_array[0]),sizeof(double)*natoms64);
   chkfile.write(reinterpret_cast<const char*>(&atoms::y_spin_array[0]),sizeof(double)*natoms64);
   chkfile.write(reinterpret_cast<const char*>(&atoms::z_spin_array[0]),sizeof(double)*natoms64);

   // Save positions and velocities only for SLD simulations
   if(sld::enabled){

      // write position array to file
      chkfile.write(reinterpret_cast<const char*>(&atoms::x_coord_array[0]),sizeof(double)*natoms64);
      chkfile.write(reinterpret_cast<const char*>(&atoms::y_coord_array[0]),sizeof(double)*natoms64);
      chkfile.write(reinterpret_cast<const char*>(&atoms::z_coord_array[0]),sizeof(double)*natoms64);

     // write velocity array to file
     chkfile.write(reinterpret_cast<const char*>(&atoms::x_velo_array[0]),sizeof(double)*natoms64);
     chkfile.write(reinterpret_cast<const char*>(&atoms::y_velo_array[0]),sizeof(double)*natoms64);
     chkfile.write(reinterpret_cast<const char*>(&atoms::z_velo_array[0]),sizeof(double)*natoms64);

  }

   // write statistical properties to file
   stats::system_magnetization.save_checkpoint(chkfile);
   stats::grain_magnetization.save_checkpoint(chkfile);
   stats::material_magnetization.save_checkpoint(chkfile);
   stats::material_grain_magnetization.save_checkpoint(chkfile);
   stats::height_magnetization.save_checkpoint(chkfile);
   stats::material_height_magnetization.save_checkpoint(chkfile);
   stats::material_grain_height_magnetization.save_checkpoint(chkfile);

   stats::system_specific_heat.save_checkpoint(chkfile);
   stats::grain_specific_heat.save_checkpoint(chkfile);
   stats::material_specific_heat.save_checkpoint(chkfile);

   stats::system_susceptibility.save_checkpoint(chkfile);
   stats::grain_susceptibility.save_checkpoint(chkfile);
   stats::material_susceptibility.save_checkpoint(chkfile);

   // write the decomposition (number of processes and region of this process); the spins are
   // stored in local atom order, so loading into another decomposition uses the atom keys written below
   const int64_t num_processors64 = vmpi::num_processors;
   chkfile.write(reinterpret_cast<const char*>(&num_processors64),sizeof(int64_t));
   chkfile.write(reinterpret_cast<const char*>(&vmpi::min_dimensions[0]),sizeof(double)*3);
   chkfile.write(reinterpret_cast<const char*>(&vmpi::max_dimensions[0]),sizeof(double)*3);

   // HAMR simulation: state of the head
   if(program::program == program::program_t::hamr_simulation) hamr::save_checkpoint(chkfile);

   // key and material of every atom, with the decomposition and the position of the spins,
   // so that the file can be loaded into another decomposition; found from the last 16 bytes
   if(create::atom_keys_valid){
      const uint64_t section = uint64_t(chkfile.tellp());
      std::vector<int32_t> materials(atoms::type_array.begin(), atoms::type_array.begin() + natoms64);
      chkfile.write(reinterpret_cast<const char*>(&num_processors64),sizeof(int64_t));
      chkfile.write(reinterpret_cast<const char*>(&vmpi::min_dimensions[0]),sizeof(double)*3);
      chkfile.write(reinterpret_cast<const char*>(&vmpi::max_dimensions[0]),sizeof(double)*3);
      chkfile.write(reinterpret_cast<const char*>(&spins_position),sizeof(uint64_t));
      chkfile.write(reinterpret_cast<const char*>(&natoms64),sizeof(uint64_t));
      chkfile.write(reinterpret_cast<const char*>(&create::atom_keys[0]),sizeof(uint64_t)*natoms64);
      chkfile.write(reinterpret_cast<const char*>(&materials[0]),sizeof(int32_t)*natoms64);
      chkfile.write(reinterpret_cast<const char*>(&section),sizeof(uint64_t));
      chkfile.write(reinterpret_cast<const char*>(&keyed_section_magic),sizeof(uint64_t));
   }

   // close checkpoint file
   chkfile.close();

   // log writing checkpoint file (only for non-continuous checkpoint files)
   if(!sim::save_checkpoint_continuous_flag) zlog << zTs() << "Checkpoint file written to disk." << std::endl;

   return;

}

//-----------------------------------------------------------------------------
// Function to save checkpoint file
//-----------------------------------------------------------------------------
void load_checkpoint(){

   // convert number of atoms, rank and time to standard long int
   uint64_t natoms64;
   int64_t time64;
   int64_t eqtime64;
   int64_t parity64;
   int64_t iH64;
   double temp;
   int64_t output_atoms_file_counter64;
   int64_t output_cells_file_counter64;
   int64_t output_rate_counter64;
   double constr_theta;
   double constr_phi;
   bool flag_constraint_theta_changed;
   bool flag_constraint_phi_changed  ;

   // variables for loading state of random number generator
   std::vector<uint32_t> mt_state(624); // 624 is hard coded in mt implementation. uint64 assumes same size as unsigned long
   int32_t mt_p=0; // position in rng state

   // determine the checkpoint file of this process
   // if there are more processes than those used when originally writing the checkpoint,
   // the data is read from the file of process 0
   std::string chkfilename = checkpoint_file_name(vmpi::my_rank);
   std::ifstream chkfile;
   chkfile.open(chkfilename.c_str(),std::ios::binary);
   bool own_file = chkfile.is_open();
   if(!own_file){
      chkfilename = checkpoint_file_name(0);
      chkfile.open(chkfilename.c_str(),std::ios::binary);
   }

   // check for open file
   if(!chkfile.is_open()){
      terminaltextcolor(RED);
      std::cerr << "Error: Unable to open checkpoint file " << chkfilename << " for reading. Exiting." << std::endl;
      std::cerr << "Info: sim:continue may be specified in the input file which requires a valid checkpoint file." << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << "Error: Unable to open checkpoint file " << chkfilename << " for reading. Exiting." << std::endl;
      zlog << zTs() << "Info: sim:continue may be specified in the input file which requires a valid checkpoint file." << std::endl;
      err::vexit();
   }

   // Set flag to true do determine that this is the beginning of the simulation
   sim::checkpoint_loaded_flag=true;
   zlog << zTs() << "Flag:checkpoint_loaded_flag = " << sim::checkpoint_loaded_flag <<std::endl;

   // read checkpoint variables from file
   chkfile.read((char*)&natoms64,sizeof(uint64_t));
   chkfile.read((char*)&time64,sizeof(int64_t));
   chkfile.read((char*)&eqtime64,sizeof(int64_t));
   chkfile.read((char*)&parity64,sizeof(int64_t));
   chkfile.read((char*)&iH64,sizeof(int64_t));
   chkfile.read((char*)&temp,sizeof(double));
   chkfile.read((char*)&constr_theta,sizeof(double));
   chkfile.read((char*)&constr_phi,sizeof(double));
   chkfile.read((char*)&flag_constraint_theta_changed,sizeof(bool));
   chkfile.read((char*)&flag_constraint_phi_changed  ,sizeof(bool));
   chkfile.read((char*)&output_atoms_file_counter64,sizeof(int64_t));
   chkfile.read((char*)&output_cells_file_counter64,sizeof(int64_t));
   chkfile.read((char*)&output_rate_counter64,sizeof(int64_t));
   chkfile.read((char*)&mt_p,sizeof(int32_t));
   chkfile.read((char*)&mt_state[0],sizeof(uint32_t)*mt_state.size());

   const uint64_t num_local_atoms = uint64_t(atoms::num_atoms-vmpi::num_halo_atoms);

   // Load saved parameters if simulation continuing
   if(sim::load_checkpoint_continue_flag){
      sim::parity = parity64;
      sim::iH = iH64;
      sim::time = time64;
      sim::equilibration_time = eqtime64;
      sim::temperature = temp;
      sim::output_atoms_file_counter = output_atoms_file_counter64;
      sim::output_cells_file_counter = output_cells_file_counter64;
      sim::output_rate_counter = output_rate_counter64;
      montecarlo::constraint_theta = constr_theta;
      montecarlo::constraint_phi = constr_phi;
      montecarlo::constraint_theta_changed = flag_constraint_theta_changed;
      montecarlo::constraint_phi_changed   = flag_constraint_phi_changed  ;
   }

   // skip the per-atom data (spins, and positions and velocities for SLD simulations),
   // which are read once the decomposition of the file is known
   const std::streampos spins_position = chkfile.tellg();
   chkfile.seekg(std::streamoff((sld::enabled ? 9 : 3)*sizeof(double)*natoms64), std::ios::cur);

   // load statistical properties from file
   stats::system_magnetization.load_checkpoint(chkfile,sim::load_checkpoint_continue_flag);
   stats::grain_magnetization.load_checkpoint(chkfile,sim::load_checkpoint_continue_flag);
   stats::material_magnetization.load_checkpoint(chkfile,sim::load_checkpoint_continue_flag);
   stats::material_grain_magnetization.load_checkpoint(chkfile,sim::load_checkpoint_continue_flag);
   stats::height_magnetization.load_checkpoint(chkfile,sim::load_checkpoint_continue_flag);
   stats::material_height_magnetization.load_checkpoint(chkfile,sim::load_checkpoint_continue_flag);
   stats::material_grain_height_magnetization.load_checkpoint(chkfile,sim::load_checkpoint_continue_flag);

   stats::system_specific_heat.load_checkpoint(chkfile,sim::load_checkpoint_continue_flag);
   stats::grain_specific_heat.load_checkpoint(chkfile,sim::load_checkpoint_continue_flag);
   stats::material_specific_heat.load_checkpoint(chkfile,sim::load_checkpoint_continue_flag);

   stats::system_susceptibility.load_checkpoint(chkfile,sim::load_checkpoint_continue_flag);
   stats::grain_susceptibility.load_checkpoint(chkfile,sim::load_checkpoint_continue_flag);
   stats::material_susceptibility.load_checkpoint(chkfile,sim::load_checkpoint_continue_flag);

   // check the decomposition (not stored in older checkpoint files, which end here)
   int64_t num_processors64 = 0;
   double min_dimensions[3] = {0.0, 0.0, 0.0};
   double max_dimensions[3] = {0.0, 0.0, 0.0};
   chkfile.read((char*)&num_processors64,sizeof(int64_t));
   chkfile.read((char*)&min_dimensions[0],sizeof(double)*3);
   chkfile.read((char*)&max_dimensions[0],sizeof(double)*3);
   const bool decomposition_stored = bool(chkfile);
   bool same = own_file && natoms64 == num_local_atoms;
   if(decomposition_stored){
      same = same && (num_processors64 == vmpi::num_processors);
      for(int i = 0; i < 3; i++) same = same && (min_dimensions[i] == vmpi::min_dimensions[i]) && (max_dimensions[i] == vmpi::max_dimensions[i]);
   }
   // all processes must take the same path (a file of this process may be left over from an earlier simulation)
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &same, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
   #endif

   // HAMR simulation: state of the head (needed only to continue)
   if(program::program == program::program_t::hamr_simulation && sim::load_checkpoint_continue_flag) hamr::load_checkpoint(chkfile);

   if(same){

      // same decomposition: the atoms are stored in local order
      chkfile.clear();
      chkfile.seekg(spins_position);
      chkfile.read((char*)&atoms::x_spin_array[0],sizeof(double)*natoms64);
      chkfile.read((char*)&atoms::y_spin_array[0],sizeof(double)*natoms64);
      chkfile.read((char*)&atoms::z_spin_array[0],sizeof(double)*natoms64);

      // Load positions and velocities only for SLD simulations
      if(sld::enabled){

         chkfile.read((char*)&atoms::x_coord_array[0],sizeof(double)*natoms64);
         chkfile.read((char*)&atoms::y_coord_array[0],sizeof(double)*natoms64);
         chkfile.read((char*)&atoms::z_coord_array[0],sizeof(double)*natoms64);

         chkfile.read((char*)&atoms::x_velo_array[0],sizeof(double)*natoms64);
         chkfile.read((char*)&atoms::y_velo_array[0],sizeof(double)*natoms64);
         chkfile.read((char*)&atoms::z_velo_array[0],sizeof(double)*natoms64);

      }

      // if continuing set state of rng
      if(sim::load_checkpoint_continue_flag) mtrandom::grnd.set_state(mt_state, mt_p);

   }
   else{

      // different decomposition: atoms are found by their keys (files written by this version)
      if(!decomposition_stored || keyed_section_position(chkfile) == 0){
         std::ostringstream message;
         if(own_file && natoms64 != num_local_atoms){
            message << "Error: Mismatch between number of atoms in checkpoint file (" << natoms64 << ") and number of generated atoms (" << num_local_atoms << "). Exiting.";
         }
         else{
            message << "Error: Checkpoint file " << chkfilename << " was written with a different decomposition (" << num_processors64
                    << " processes) and holds no atom keys. Restart with the same number of processes as the simulation that wrote it. Exiting.";
         }
         err::fatal(message.str());
      }
      if(sld::enabled) err::fatal("Error: Spin-lattice simulations can only be restarted with the same decomposition as the simulation that wrote the checkpoint. Exiting.");
      if(!create::atom_keys_valid) err::fatal("Error: The atoms cannot be identified by their unit cell, which loading a checkpoint written with a different decomposition requires. Exiting.");

      load_keyed_spins();
      zlog << zTs() << "Checkpoint loaded by atom key." << std::endl;

      // the random number states of the processes that wrote the checkpoint do not map onto
      // these processes: continue with new states, derived from the seed and the time step
      if(sim::load_checkpoint_continue_flag){
         mtrandom::grnd.seed(vmpi::parallel_rng_seed(int(uint32_t(mtrandom::integration_seed) + 1000003u * uint32_t(sim::time))));
         zlog << zTs() << "CPU random number generator reseeded at time step " << sim::time << ": its random numbers differ from those of an uninterrupted simulation." << std::endl;
      }

   }

   // the file holds the spins of local atoms only: update the halo atoms, which
   // statistics calculated before the first time step use
   #ifdef MPICF
      vmpi::mpi_init_halo_swap();
      vmpi::mpi_complete_halo_swap();
   #endif

   // close checkpoint file
   chkfile.close();

   // log reading checkpoint file
   zlog << zTs() << "Checkpoint file loaded at sim::time " << sim::time << "." << std::endl;

   return;

}
