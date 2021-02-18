//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2021. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
// Vampire headers
#include "micromagnetic.hpp"
#include "sim.hpp"

// micromagnetic module headers
#include "internal.hpp"
//#include "material.hpp"

// shorthand for brevity
namespace mm = micromagnetic::internal;

//#include "cells.hpp"

namespace micromagnetic{
   namespace internal{

      //-----------------------------------------------------------------------------
      // calculates the stt precession term for each cell from the atomic parameters
      //-----------------------------------------------------------------------------
      void calculate_stt(const int num_atoms,
                         const int num_cells,
                         const std::vector<int>& cell_array,      // 1D array storing which cell each atom is in
                         const std::vector<int>& type_array,      // 1D array storing which material each atom is
                         std::vector <mp::materials_t>& material, //class of material parameters for the atoms
                         std::vector<double>& stt_rj,
                         std::vector<double>& stt_pj){

         // resize arrays storing cell level properties
         stt_pj.resize(num_cells,0.0);
         stt_rj.resize(num_cells,0.0);

         //vector stores the number of atoms per micromagnetic cell
         std::vector<double> atoms_pc(num_cells,0.0);                                        //vector stores the number of atoms per micromagnetic cell

         // get material level stt constants from sim module
         std::vector<double> stt_rj_mat = sim::get_stt_rj();
         std::vector<double> stt_pj_mat = sim::get_stt_pj();

         //sums over all local atoms to add the stt vaues per cell
         for (int atom = 0; atom < num_atoms; atom++){
            int cell = cell_array[atom];
            int mat = type_array[atom];
            stt_rj[cell] += stt_rj_mat[mat];
            stt_pj[cell] += stt_pj_mat[mat];
            atoms_pc[cell]++;
         }

         // In parallel reduce stt parameters on all cells
         #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, &stt_rj[0],   num_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &stt_pj[0],   num_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &atoms_pc[0], num_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         #endif

         // normalise the total stt parameters on all cells (on all processors)
         for (int cell = 0; cell < num_cells; cell++){
            stt_rj[cell] = stt_rj[cell] / atoms_pc[cell];
            stt_pj[cell] = stt_pj[cell] / atoms_pc[cell];
         }

         // initialise spin transfer torque polarization
         std::vector<double> sttp = sim::get_stt_polarization_unit_vector();

         mm::sttpx = sttp[0];
         mm::sttpy = sttp[1];
         mm::sttpz = sttp[2];

         return;

      }

   } //closes the internal namspace
}  //closes the micromagnetic namespace
