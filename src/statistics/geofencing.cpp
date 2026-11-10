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
//   Geofencing identifies the magnetic state of a magnetisation direction m as
//   the closest easy axis e_a, if m . e_a >= threshold. Otherwise m is lost
//   between states. The easy axes are the eight <111> or six <100> directions
//   for cubic geofencing, or the positive and negative easy axis for uniaxial
//   geofencing. States are numbered from 1 in the order of the axes below.
//
//   The state of the system or material magnetisation is output directly. For
//   per-spin geofencing every spin is tracked separately for each threshold,
//   and the mean over spins is output of the switching time
//
//      tau = (time of last transition - time of first transition) / transitions
//
//   (zero for a spin with no transitions), the percentage of samples in which
//   a spin is lost, and the total number of transitions.
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>

// Vampire headers
#include "errors.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// statistics module headers
#include "internal.hpp"

namespace stats{

   namespace internal{

      //------------------------------------------------------------------------
      // Function to find the closest easy axis to the direction (mx, my, mz),
      // returning the axis index and m . e_a for that axis
      //------------------------------------------------------------------------
      void closest_axis(const std::vector<double>& axes, const double mx, const double my, const double mz, int& best_axis, double& best_dot){

         best_axis = 0;
         best_dot = mx*axes[0] + my*axes[1] + mz*axes[2];

         const int num_axes = axes.size()/3;
         for(int axis = 1; axis < num_axes; axis++){
            const double dot = mx*axes[3*axis] + my*axes[3*axis+1] + mz*axes[3*axis+2];
            if(dot > best_dot){
               best_dot = dot;
               best_axis = axis;
            }
         }

         return;

      }

      //------------------------------------------------------------------------
      // Function to set the easy axes, check the materials and set up per-spin
      // geofencing for the local spins included in the statistics
      //------------------------------------------------------------------------
      void initialize_geofencing(const int num_atoms, const int num_materials, const std::vector<bool>& magnetic){

         // cubic easy axes, <111> or <100> when the <111> directions are hard
         const double a = 1.0/sqrt(3.0);
         if(cubic_geofencing_111_hard) geofencing[cubic_geofencing].axes = { 1,0,0, -1,0,0, 0,1,0, 0,-1,0, 0,0,1, 0,0,-1 };
         else geofencing[cubic_geofencing].axes = { a,a,a, -a,a,a, a,-a,a, -a,-a,a, a,a,-a, -a,a,-a, a,-a,-a, -a,-a,-a };

         // uniaxial easy axes, parallel and antiparallel to the easy axis
         const std::vector<double>& u = uniaxial_geofencing_axis;
         geofencing[uniaxial_geofencing].axes = { u[0], u[1], u[2], -u[0], -u[1], -u[2] };

         // check that geofenced materials exist
         if(geofencing[cubic_geofencing].material >= num_materials) err::fatal("Error - output:geofencing-cubic-material-id is greater than the number of materials (" + std::to_string(num_materials) + ")");
         if(geofencing[uniaxial_geofencing].material >= num_materials) err::fatal("Error - output:geofencing-uniaxial-material-id is greater than the number of materials (" + std::to_string(num_materials) + ")");

         if(!calculate_cubic_per_spin_geofencing && !calculate_uniaxial_per_spin_geofencing) return;

         // include only magnetic spins and count them on all processors
         per_spin_geofencing_mask.assign(magnetic.begin(), magnetic.begin() + num_atoms);
         num_per_spin_geofencing_atoms = std::count(per_spin_geofencing_mask.begin(), per_spin_geofencing_mask.end(), true);
         #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, &num_per_spin_geofencing_atoms, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
         #endif

         if(calculate_cubic_per_spin_geofencing) reset_per_spin_geofencing(cubic_geofencing);
         if(calculate_uniaxial_per_spin_geofencing) reset_per_spin_geofencing(uniaxial_geofencing);

         return;

      }

      //------------------------------------------------------------------------
      // Function to geofence each spin and record its transitions between states
      //------------------------------------------------------------------------
      void update_per_spin_geofencing(const geofencing_t type, const std::vector<double>& sx, const std::vector<double>& sy, const std::vector<double>& sz){

         geofencing_data_t& g = geofencing[type];
         const int num_thresholds = per_spin_geofencing_thresholds.size();
         const int num_atoms = per_spin_geofencing_mask.size();
         const double time = sim::time*mp::dt_SI; // current time (s)

         g.num_samples++;

         for(int atom = 0; atom < num_atoms; atom++){

            if(!per_spin_geofencing_mask[atom]) continue;

            int axis;
            double dot;
            closest_axis(g.axes, sx[atom], sy[atom], sz[atom], axis, dot);

            for(int t = 0; t < num_thresholds; t++){

               // spin is lost between states, sums are stored as [tau, lost steps, transitions] for each threshold
               if(dot < per_spin_geofencing_thresholds[t]){
                  g.sums[num_thresholds + t] += 1.0;
                  continue;
               }

               // transition from the previous state
               int& state = g.state[t][atom];
               if(state >= 0 && axis != state){

                  uint64_t& transitions = g.transitions[t][atom];
                  double& first_time = g.first_transition_time[t][atom];
                  double& last_time = g.last_transition_time[t][atom];

                  // replace tau for this spin in the sum of tau over spins
                  if(transitions > 0) g.sums[t] -= (last_time - first_time)/double(transitions);
                  else first_time = time;
                  last_time = time;
                  transitions++;
                  g.sums[t] += (last_time - first_time)/double(transitions);

                  g.sums[2*num_thresholds + t] += 1.0;

               }
               state = axis;

            }

         }

         // sum over all processors
         g.totals = g.sums;
         #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, g.totals.data(), g.totals.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         #endif

         return;

      }

      //------------------------------------------------------------------------
      // Function to reset per-spin geofencing data
      //------------------------------------------------------------------------
      void reset_per_spin_geofencing(const geofencing_t type){

         geofencing_data_t& g = geofencing[type];
         const int num_thresholds = per_spin_geofencing_thresholds.size();
         const int num_atoms = per_spin_geofencing_mask.size();

         g.state.assign(num_thresholds, std::vector<int>(num_atoms, -1));
         g.transitions.assign(num_thresholds, std::vector<uint64_t>(num_atoms, 0));
         g.first_transition_time.assign(num_thresholds, std::vector<double>(num_atoms, 0.0));
         g.last_transition_time.assign(num_thresholds, std::vector<double>(num_atoms, 0.0));
         g.sums.assign(3*num_thresholds, 0.0);
         g.totals.assign(3*num_thresholds, 0.0);
         g.num_samples = 0;

         return;

      }

      //------------------------------------------------------------------------
      // Function to format a threshold for column headers, e.g. 0.950
      //------------------------------------------------------------------------
      std::string threshold_label(const double threshold){
         std::ostringstream label;
         label.setf(std::ios::fixed, std::ios::floatfield);
         label.precision(3);
         label << threshold;
         return label.str();
      }

   } // end of internal namespace

   //---------------------------------------------------------------------------
   // Function to output the geofencing state of the system or material
   // magnetisation, and optionally m . e_a for the closest axis. A lost state
   // is output as 0X, where X is the closest state.
   //---------------------------------------------------------------------------
   std::string output_geofencing_state(const geofencing_t type, const bool header){

      const internal::geofencing_data_t& g = internal::geofencing[type];
      const bool output_product = (type == cubic_geofencing) ? output_cubic_geofencing_product : output_uniaxial_geofencing_product;
      const std::string name = std::string(type == cubic_geofencing ? "cubic" : "uniaxial") + "_geofencing_";
      const std::string source = (g.material < 0) ? "global" : "m" + std::to_string(g.material + 1);

      // result string stream
      std::ostringstream res;

      // set custom precision if enabled
      if(vout::custom_precision){
         res.precision(vout::precision);
         if(vout::fixed) res.setf(std::ios::fixed, std::ios::floatfield);
      }
      vout::fixed_width_output state_column(res, vout::fw_size_int);
      vout::fixed_width_output product_column(res, vout::fw_size);

      if(header){
         state_column << name + source;
         if(output_product) product_column << name + "dot_" + source;
         return res.str();
      }

      // magnetisation direction [mx, my, mz, |m|] for the system or material
      const std::vector<double>& m = (g.material < 0) ? system_magnetization.get_magnetization() : material_magnetization.get_magnetization();
      const int index = (g.material < 0) ? 0 : 4*g.material;

      int axis;
      double dot;
      internal::closest_axis(g.axes, m[index], m[index+1], m[index+2], axis, dot);

      if(dot >= g.threshold) state_column << axis + 1;
      else state_column << "0" + std::to_string(axis + 1);
      if(output_product) product_column << dot;

      return res.str();

   }

   //---------------------------------------------------------------------------
   // Function to output per-spin geofencing data for each threshold: mean tau (s),
   // mean percentage of samples lost and total number of transitions
   //---------------------------------------------------------------------------
   std::string output_per_spin_geofencing(const geofencing_t type, const bool header){

      const internal::geofencing_data_t& g = internal::geofencing[type];
      const std::vector<double>& thresholds = internal::per_spin_geofencing_thresholds;
      const int num_thresholds = thresholds.size();
      const std::string name = std::string("PerSpin_") + (type == cubic_geofencing ? "Cubic" : "Uniaxial") + "_Geofencing_";

      // result string stream
      std::ostringstream res;

      // set custom precision if enabled
      if(vout::custom_precision){
         res.precision(vout::precision);
         if(vout::fixed) res.setf(std::ios::fixed, std::ios::floatfield);
      }
      vout::fixed_width_output double_column(res, vout::fw_size);
      vout::fixed_width_output int_column(res, vout::fw_size_int);

      if(header){
         for(int t = 0; t < num_thresholds; t++) double_column << name + "Tau_" + internal::threshold_label(thresholds[t]);
         for(int t = 0; t < num_thresholds; t++) double_column << name + "LostTime_" + internal::threshold_label(thresholds[t]);
         for(int t = 0; t < num_thresholds; t++) int_column << name + "Transitions_" + internal::threshold_label(thresholds[t]);
         return res.str();
      }

      // means over all spins
      const double num_atoms = internal::num_per_spin_geofencing_atoms;
      const double inv_atoms = (num_atoms > 0) ? 1.0/num_atoms : 0.0;
      const double inv_samples = (g.num_samples > 0) ? 1.0/double(g.num_samples) : 0.0;

      for(int t = 0; t < num_thresholds; t++) double_column << g.totals[t]*inv_atoms;
      for(int t = 0; t < num_thresholds; t++) double_column << 100.0*g.totals[num_thresholds + t]*inv_atoms*inv_samples;
      for(int t = 0; t < num_thresholds; t++) int_column << uint64_t(g.totals[2*num_thresholds + t]);

      return res.str();

   }

} // end of stats namespace
