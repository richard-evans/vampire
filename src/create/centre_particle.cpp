//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "errors.hpp"
#include "create.hpp"
#include "vmpi.hpp"

// Internal create header
#include "internal.hpp"

namespace create{
namespace internal{

//------------------------------------------------------------------------------------------------------
// Function to alter particle origin to be centred on an atom
//------------------------------------------------------------------------------------------------------
void centre_particle_on_atom(std::vector<double>& particle_origin, std::vector<cs::catom_t>& catom_array){

   vmpi::barrier();

   // set initial max range
   double max_range_sq = 1e123;
   unsigned int nearest; // nearest atom to initial particle origin

   // copy to temporary for speed
   const double prx = particle_origin[0];
   const double pry = particle_origin[1];
   const double prz = particle_origin[2];

   // loop over all atoms to find closest atom
   for(size_t atom=0;atom<catom_array.size();atom++){
      double dx = catom_array[atom].x - prx;
      double dy = catom_array[atom].y - pry;
      double dz = catom_array[atom].z - prz;
      double r = dx*dx + dy*dy + dz*dz;
      if(r < max_range_sq){
         max_range_sq = r;
         nearest = atom;
      }
   }

   // set particle origin to nearest atom
   particle_origin[0] = catom_array[nearest].x;
   particle_origin[1] = catom_array[nearest].y;
   particle_origin[2] = catom_array[nearest].z;

   //-----------------------------------------------------
   // For parallel reduce on all CPUs
   //-----------------------------------------------------
   #ifdef MPICF

      // set up array to get ranges from all CPUs on rank 0
      std::vector<double> ranges;
      if(vmpi::my_rank == 0) ranges.resize(vmpi::num_processors, 1.e123);
      else ranges.resize(1,0.0); // one value sufficient on all other CPUs

      // gather max ranges from all cpus on root (1 data point from each process)
      MPI_Gather(&max_range_sq, 1, MPI_DOUBLE, &ranges[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      // variable to store rank of minimum range
      unsigned int rank_of_min_range=0;

      // work out minimum range on root
      if(vmpi::my_rank==0){

         double min_range = 1.e123;

         // loop over all ranges and determine minimum and cpu location
         for(int i=0; i<ranges.size(); i++){
            if(ranges[i] < min_range){
               min_range = ranges[i];
               rank_of_min_range = i;
            }
         }
      }

      // broadcast id of nearest to all cpus from root
      MPI_Bcast(&rank_of_min_range, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

      // broadcast position to all cpus
      MPI_Bcast(&particle_origin[0], 3, MPI_DOUBLE, rank_of_min_range, MPI_COMM_WORLD);

      vmpi::barrier();

   #endif

   return;

}

} // end of namespace internal
} // end of namespace create
