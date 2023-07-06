//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// System headers
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

// Program headers
#include "atoms.hpp"
#include "errors.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "program.hpp"

//-----------------------------------------------------------------------------
// Function to save checkpoint file
//-----------------------------------------------------------------------------
void save_checkpoint(){

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
   double constr_theta = sim::constraint_theta;
   double constr_phi   = sim::constraint_phi;
   bool flag_constraint_theta_changed = sim::constraint_theta_changed;
   bool flag_constraint_phi_changed   = sim::constraint_phi_changed;

   // determine checkpoint file name
   std::stringstream chkfilenamess;
   chkfilenamess << "vampire" << vmpi::my_rank << ".chk";
   std::string chkfilename = chkfilenamess.str();

   // open checkpoint file
   std::ofstream chkfile;
   chkfile.open(chkfilename.c_str(),std::ios::binary);

   // check for open file
   if(!chkfile.is_open()){
      terminaltextcolor(RED);
      std::cerr << "Error: Unable to open checkpoint file " << chkfilename << " for writing. Exiting." << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << "Error: Unable to open checkpoint file " << chkfilename << " for writing. Exiting." << std::endl;
      err::vexit();
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
   chkfile.write(reinterpret_cast<const char*>(&atoms::x_spin_array[0]),sizeof(double)*natoms64);
   chkfile.write(reinterpret_cast<const char*>(&atoms::y_spin_array[0]),sizeof(double)*natoms64);
   chkfile.write(reinterpret_cast<const char*>(&atoms::z_spin_array[0]),sizeof(double)*natoms64);

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

   // determine checkpoint file name
   std::stringstream chkfilenamess;
   chkfilenamess << "vampire" << vmpi::my_rank << ".chk";
   std::string chkfilename = chkfilenamess.str();

   // open checkpoint file
   std::ifstream chkfile;
   chkfile.open(chkfilename.c_str(),std::ios::binary);

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

   //std::cout << "random generator state loaded = " << mt_p << std::endl;
   // if continuing set state of rng
   if(sim::load_checkpoint_continue_flag) mtrandom::grnd.set_state(mt_state, mt_p);

   // check for rational number of atoms
   if(static_cast<uint64_t>(atoms::num_atoms-vmpi::num_halo_atoms) != natoms64){
      terminaltextcolor(RED);
      std::cerr << "Error: Mismatch between number of atoms in checkpoint file (" << natoms64 << ") and number of generated atoms (" << atoms::num_atoms-vmpi::num_halo_atoms << "). Exiting." << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << "Error: Mismatch between number of atoms in checkpoint file (" << natoms64 << ") and number of generated atoms (" << atoms::num_atoms-vmpi::num_halo_atoms << "). Exiting." << std::endl;
      err::vexit();
   }

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
      sim::constraint_theta = constr_theta;
      sim::constraint_phi = constr_phi;
      sim::constraint_theta_changed = flag_constraint_theta_changed;
      sim::constraint_phi_changed   = flag_constraint_phi_changed  ;
   }

   // Load spin positions
   chkfile.read((char*)&atoms::x_spin_array[0],sizeof(double)*natoms64);
   chkfile.read((char*)&atoms::y_spin_array[0],sizeof(double)*natoms64);
   chkfile.read((char*)&atoms::z_spin_array[0],sizeof(double)*natoms64);

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

   // close checkpoint file
   chkfile.close();

   // log reading checkpoint file
   zlog << zTs() << "Checkpoint file loaded at sim::time " << sim::time << "." << std::endl;

   return;

}
