//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Lennon Sivaruban 2026. All rights reserved.
//
//   Email: lennon.sivaruban@manchester.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef STATS_INTERNAL_H_
#define STATS_INTERNAL_H_

// C++ standard library headers
#include <cstdint>
#include <vector>

// Vampire headers
#include "stats.hpp"

namespace stats{

   namespace internal{

      //------------------------------------------------------------------------
      // Settings and per-spin data for one type of geofencing (cubic or uniaxial)
      //------------------------------------------------------------------------
      struct geofencing_data_t{

         std::vector<double> axes; // easy axes e_a, 3 components per axis (state = axis index + 1)
         double threshold = 0.95;  // magnetisation is in a state if m . e_a >= threshold
         int material = -1;        // material magnetisation to geofence (-1 for system magnetisation)

         // per-spin data for each threshold and local atom [threshold][atom]
         std::vector<std::vector<int>> state;                     // current state of each spin (-1 before the first state)
         std::vector<std::vector<uint64_t>> transitions;          // number of transitions of each spin between states
         std::vector<std::vector<double>> first_transition_time;  // time of first transition (s)
         std::vector<std::vector<double>> last_transition_time;   // time of last transition (s)

         // per-spin sums over local spins for each threshold [tau, lost steps, transitions], and their totals over all processors
         std::vector<double> sums;
         std::vector<double> totals;
         uint64_t num_samples = 0; // number of samples since last reset

      };

      extern geofencing_data_t geofencing[2]; // geofencing data indexed by stats::geofencing_t

      extern bool cubic_geofencing_111_hard;                     // use six <100> easy axes for cubic geofencing instead of eight <111> axes
      extern std::vector<double> uniaxial_geofencing_axis;       // easy axis for uniaxial geofencing
      extern std::vector<double> per_spin_geofencing_thresholds; // m . e_a thresholds for per-spin geofencing
      extern std::vector<bool> per_spin_geofencing_mask;         // spins included in per-spin geofencing (magnetic atoms)
      extern uint64_t num_per_spin_geofencing_atoms;             // total number of spins in per-spin geofencing

      //------------------------------------------------------------------------
      // Internal function declarations
      //------------------------------------------------------------------------
      void initialize_geofencing(const int num_atoms, const int num_materials, const std::vector<bool>& magnetic);
      void update_per_spin_geofencing(const geofencing_t type, const std::vector<double>& sx, const std::vector<double>& sy, const std::vector<double>& sz);
      void reset_per_spin_geofencing(const geofencing_t type);

   } // end of internal namespace

} // end of stats namespace

#endif //STATS_INTERNAL_H_
