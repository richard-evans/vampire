//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   GNU GPL (version 2) licence (see licence file for details).
//
//   (c) S Jenkins 2026. All rights reserved.
//
//------------------------------------------------------------------------------
//
//   Active window for the HAMR write (hamr:active-window).
//
//   Without a dipole field, and with grains that share no exchange bonds, each
//   grain evolves only from its own spins, its laser temperature and the head
//   field. Far from the head a grain sits at the minimum temperature with no
//   field, where it cannot switch on the time scale of the write, so it is
//   frozen: its spins are kept and not integrated. A grain is integrated while
//   the laser raises its temperature anywhere by at least
//   hamr:active-window-temperature above the minimum, or it overlaps the head
//   field box, and for a further relaxation time (10 ps) after that, so that it
//   is frozen in equilibrium at the minimum temperature (also when the head
//   jumps to the next track). All grains are integrated for the first 10 ps.
//   The window is not changed while the laser is off.
//
//   Window snapshots: while the window runs, only the atoms of grains
//   integrated since the previous snapshot can have changed, so a spin
//   snapshot may hold just those atoms (src/config/atoms_spins.cpp).
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

// Vampire headers
#include "atoms.hpp"
#include "anisotropy.hpp"
#include "dipole.hpp"
#include "errors.hpp"
#include "exchange.hpp"
#include "gpu.hpp"
#include "hamr.hpp"
#include "integrate.hpp"
#include "material.hpp"
#include "program.hpp"
#include "sim.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// hamr headers
#include "internal.hpp"

namespace hamr{

   namespace internal{

      // per grain: bounding box of its magnetic atoms, steps since it was last hot, integrated or not
      std::vector<double> grain_min_x, grain_max_x, grain_min_y, grain_max_y;
      std::vector<uint64_t> grain_cool_steps;
      std::vector<bool> grain_active;
      std::vector<int> active_ranges; // [first, last+1) pairs of integrated atoms
      std::vector<char> grain_changed;  // grain integrated since the last spin snapshot
      bool snapshot_in_window = false;  // a spin snapshot was written while the window was running
      uint64_t window_hold_steps = 0;
      uint64_t num_window_changes = 0;
      bool window_initialised = false;
      std::vector<uint64_t> resume_grain_cool_steps; // steps since each grain was last hot, from a checkpoint

      // marks the window state in a checkpoint file (older files end the HAMR state before it)
      const uint64_t window_state_magic = 0x3157444e49575248ULL;

      bool is_magnetic(const int atom){
         return mp::material[atoms::type_array[atom]].non_magnetic == 0;
      }

      // atoms integrated by this process (halo atoms are copies of other processes' atoms)
      int num_owned_atoms(){
         #ifdef MPICF
            return vmpi::num_core_atoms + vmpi::num_bdry_atoms;
         #else
            return atoms::num_atoms;
         #endif
      }

      void rebuild_active_ranges(){
         num_window_changes++;
         active_ranges.clear();
         bool in_range = false;
         const int num_local = num_owned_atoms();
         for(int atom = 0; atom < num_local; atom++){
            const bool active = is_magnetic(atom) && grain_active[atoms::grain_array[atom]];
            if(active && !in_range) active_ranges.push_back(atom);
            if(!active && in_range) active_ranges.push_back(atom);
            in_range = active;
         }
         if(in_range) active_ranges.push_back(num_local);

         int num_grains = 0;
         for(size_t g = 0; g < grain_active.size(); g++) if(grain_active[g]) num_grains++;
         // every process changes the window at the same step, so the total can be summed here
         uint64_t num_atoms = 0;
         for(size_t r = 0; r < active_ranges.size(); r += 2) num_atoms += active_ranges[r+1] - active_ranges[r];
         num_atoms = vmpi::all_reduce_sum(num_atoms);
         zlog << zTs() << "HAMR active window: " << num_grains << " grains, " << num_atoms << " atoms integrated, head at "
              << head_position_x*0.1 << " nm, " << head_position_y*0.1 << " nm" << std::endl;
      }

      //------------------------------------------------------------------------
      // Refuses simulations the physics argument above does not cover and sets
      // up the grain bounding boxes; called when the HAMR write starts
      //------------------------------------------------------------------------
      void initialize_active_window(){

         std::vector<std::string> unsupported;
         if(integrate::integrator != integrate::integrator_t::llg_heun) unsupported.push_back("integrators other than llg-heun");
         if(::dipole::activated) unsupported.push_back("dipole fields");
         if(::exchange::get_exchange_type() != 0) unsupported.push_back("vectorial or tensorial exchange");
         if(::exchange::biquadratic || ::exchange::four_spin) unsupported.push_back("biquadratic or four spin exchange");
         if(anisotropy::is_neel_enabled()) unsupported.push_back("Neel anisotropy (couples neighbouring atoms)");
         if(program::lagrange_multiplier_enabled) unsupported.push_back("the Lagrange multiplier");
         if(sim::spin_torque_fields_enabled() || sim::vcma_fields_enabled()) unsupported.push_back("spin torque or voltage controlled anisotropy fields");

         // grains must evolve independently: no exchange between a magnetic atom and an atom outside its grain
         // (same tolerance as the GPU neighbour list); the neighbours of local atoms include halo atoms
         const double tol = gpu::exchange_tolerance;
         const int num_local = num_owned_atoms();
         int bond_i = -1, bond_j = -1;
         for(int atom = 0; atom < num_local && bond_i < 0; atom++){
            if(!is_magnetic(atom)) continue;
            for(int nn = atoms::neighbour_list_start_index[atom]; nn <= atoms::neighbour_list_end_index[atom]; nn++){
               const int j = atoms::neighbour_list_array[nn];
               if(std::fabs(atoms::i_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij) <= tol) continue;
               if(!is_magnetic(j) || atoms::grain_array[j] != atoms::grain_array[atom]){ bond_i = atom; bond_j = j; break; }
            }
         }
         int bond_found = bond_i >= 0;
         #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, &bond_found, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
         #endif
         if(bond_found){
            std::ostringstream bond;
            bond << "exchange between atoms of different grains or with non-magnetic atoms";
            if(bond_i >= 0) bond << " (atoms " << bond_i << " and " << bond_j << ")";
            unsupported.push_back(bond.str());
         }

         if(!unsupported.empty()){
            std::string message = "Error: hamr:active-window does not support the following features of this simulation:\n";
            for(size_t i = 0; i < unsupported.size(); i++) message += "\t" + unsupported[i] + "\n";
            message += "Run without hamr:active-window to use them.";
            err::fatal(message);
         }

         // grains and their boxes over all processes, so that every process makes the same decisions
         int num_grains = 0;
         for(int atom = 0; atom < num_local; atom++) if(is_magnetic(atom)) num_grains = std::max(num_grains, atoms::grain_array[atom] + 1);
         #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, &num_grains, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
         #endif
         grain_min_x.assign(num_grains,  1e300); grain_max_x.assign(num_grains, -1e300);
         grain_min_y.assign(num_grains,  1e300); grain_max_y.assign(num_grains, -1e300);
         // shortcut: a grain split by a periodic boundary gets a box across the film and is always integrated
         for(int atom = 0; atom < num_local; atom++){
            if(!is_magnetic(atom)) continue;
            const int g = atoms::grain_array[atom];
            grain_min_x[g] = std::min(grain_min_x[g], atom_coords_x[atom]);
            grain_max_x[g] = std::max(grain_max_x[g], atom_coords_x[atom]);
            grain_min_y[g] = std::min(grain_min_y[g], atom_coords_y[atom]);
            grain_max_y[g] = std::max(grain_max_y[g], atom_coords_y[atom]);
         }
         #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, grain_min_x.data(), num_grains, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, grain_min_y.data(), num_grains, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, grain_max_x.data(), num_grains, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, grain_max_y.data(), num_grains, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
         #endif

         window_hold_steps = uint64_t(round(1.0e-11/mp::dt_SI));
         grain_cool_steps.assign(num_grains, 0);
         grain_active.assign(num_grains, true);
         grain_changed.assign(num_grains, 1);
         snapshot_in_window = false;
         window_initialised = true;

         // continued simulation: the grains integrated at the checkpoint, as in the simulation that wrote it
         if(!resume_grain_cool_steps.empty()){
            if(resume_grain_cool_steps.size() == grain_cool_steps.size()){
               grain_cool_steps = resume_grain_cool_steps;
               for(int g = 0; g < num_grains; g++) grain_active[g] = grain_cool_steps[g] < window_hold_steps;
               zlog << zTs() << "HAMR active window state loaded from the checkpoint" << std::endl;
            }
            else{
               err::warning("Warning - the checkpoint holds the active window of " + std::to_string(resume_grain_cool_steps.size()) +
                            " grains but the system has " + std::to_string(num_grains) + ", so every grain is integrated for the first 10 ps");
            }
         }

         zlog << zTs() << "HAMR active window enabled: " << num_grains << " grains, temperature threshold " << active_window_dT
              << " K above the minimum, grains kept for " << window_hold_steps << " steps after leaving it" << std::endl;
         rebuild_active_ranges();
      }

      //------------------------------------------------------------------------
      // Functions to write and read the window state (steps since each grain was
      // last hot) after the head state in a checkpoint file, so that a continued
      // simulation integrates the same grains as an uninterrupted one
      //------------------------------------------------------------------------
      void save_window_state(std::ostream& chkfile){
         if(!window_initialised) return;
         const uint64_t num_grains = grain_cool_steps.size();
         chkfile.write(reinterpret_cast<const char*>(&window_state_magic), sizeof(uint64_t));
         chkfile.write(reinterpret_cast<const char*>(&num_grains), sizeof(uint64_t));
         chkfile.write(reinterpret_cast<const char*>(grain_cool_steps.data()), sizeof(uint64_t)*num_grains);
      }

      void load_window_state(std::istream& chkfile){
         uint64_t magic = 0;
         uint64_t num_grains = 0;
         chkfile.read(reinterpret_cast<char*>(&magic), sizeof(uint64_t));
         chkfile.read(reinterpret_cast<char*>(&num_grains), sizeof(uint64_t));
         if(!chkfile || magic != window_state_magic) return;
         std::vector<uint64_t> cool_steps(num_grains);
         chkfile.read(reinterpret_cast<char*>(cool_steps.data()), sizeof(uint64_t)*num_grains);
         if(chkfile) resume_grain_cool_steps.swap(cool_steps);
      }

   } // end of internal namespace

   bool active_window_requested(){
      return internal::active_window && program::program == program::program_t::hamr_simulation;
   }

   bool active_window_enabled(){
      return internal::window_initialised;
   }

   //---------------------------------------------------------------------------
   // Spin snapshots with the window: a snapshot may hold only the atoms whose
   // spins changed since the previous snapshot if that snapshot was written
   // while the window was running (every atom moves before the window starts)
   //---------------------------------------------------------------------------
   bool window_snapshot_possible(){
      return internal::window_initialised && internal::snapshot_in_window && !internal::active_window_full_snapshots;
   }

   bool spin_changed_since_snapshot(const int atom){
      return internal::is_magnetic(atom) && internal::grain_changed[atoms::grain_array[atom]];
   }

   void snapshot_written(){
      if(!internal::window_initialised) return;
      internal::snapshot_in_window = true;
      std::fill(internal::grain_changed.begin(), internal::grain_changed.end(), 0);
   }

   uint64_t window_changes(){
      return internal::num_window_changes;
   }

   bool atom_in_window(const int atom){
      return internal::is_magnetic(atom) && internal::grain_active[atoms::grain_array[atom]];
   }

   // marks the grains integrated at the current step as changed since the last snapshot
   void mark_changed_grains(){
      for(size_t g = 0; g < internal::grain_active.size(); g++) if(internal::grain_active[g]) internal::grain_changed[g] = 1;
   }

   //---------------------------------------------------------------------------
   // Updates the window for the current head position (while the laser is on)
   // and returns the [first, last+1) pairs of the atoms to integrate
   //---------------------------------------------------------------------------
   const std::vector<int>& update_active_window(){

      using namespace hamr::internal;
      if(!head_laser_on){
         mark_changed_grains();
         return active_ranges;
      }

      const double px = head_position_x;
      const double py = head_position_y;
      const double one_over_denx = 1.0/(2.0*laser_sigma_x*laser_sigma_x);
      const double one_over_deny = 1.0/(2.0*laser_sigma_y*laser_sigma_y);
      const double DeltaT = Tmax - Tmin;
      // head field box, as in apply_field_spatial_box
      const double box_min_x = px - 0.5*H_bounds_x - NPS;
      const double box_max_x = px + 0.5*H_bounds_x - NPS;
      const double box_min_y = py - 0.5*H_bounds_y;
      const double box_max_y = py + 0.5*H_bounds_y;

      bool changed = false;
      for(size_t g = 0; g < grain_active.size(); g++){
         // the laser temperature is largest at the point of the bounding box nearest the head
         const double dx = std::max(0.0, std::max(grain_min_x[g] - px, px - grain_max_x[g]));
         const double dy = std::max(0.0, std::max(grain_min_y[g] - py, py - grain_max_y[g]));
         const bool hot = DeltaT * exp(-dx*dx*one_over_denx) * exp(-dy*dy*one_over_deny) >= active_window_dT;
         const bool in_field = grain_max_x[g] >= box_min_x && grain_min_x[g] <= box_max_x &&
                               grain_max_y[g] >= box_min_y && grain_min_y[g] <= box_max_y;
         if(hot || in_field) grain_cool_steps[g] = 0;
         else if(grain_cool_steps[g] < window_hold_steps) grain_cool_steps[g]++;
         const bool active = grain_cool_steps[g] < window_hold_steps;
         if(active != grain_active[g]){ grain_active[g] = active; changed = true; }
      }
      if(changed) rebuild_active_ranges();
      mark_changed_grains();

      return active_ranges;
   }

} // end of hamr namespace
