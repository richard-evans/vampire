//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <iostream>

// Vampire headers
#include "ltmp.hpp"
#include "vmpi.hpp"

// Local temperature pulse headers
#include "internal.hpp"

namespace ltmp{
   namespace internal{

      //-----------------------------------------------------------------------------
      // Function to calculate the local temperature using the two temperature model
      //-----------------------------------------------------------------------------
      void calculate_local_temperature_pulse(const double time_from_start){

         const double reduced_time =  (time_from_start-2.*ltmp::internal::pump_time)/(ltmp::internal::pump_time);
         const double prefactor = 4.0*log(2.0); // normalise to unit width
         const double pump=ltmp::internal::pump_power*exp(-reduced_time*reduced_time*prefactor);

         const double G  = ltmp::internal::TTG;
         const double Ce = ltmp::internal::TTCe;
         const double Cl = ltmp::internal::TTCl;
         const double dt = ltmp::internal::dt;

         // Parallisation
         // if vertical only
         // loop over internal cells
         // broadcast result
         // else
         // loop over 1/n cells
         // calculate dTe dTp from Te, Tp
         //#ifdef MPICF
         //   MPI_Allreduce(MPI_IN_PLACE, &st::internal::spin_torque[0],st::internal::spin_torque.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         //#endif

         // Precalculate heat transfer constant k*L/V (J/K/m^3/s) (divide by Angstroms^2)
         const double dTdiff_prefactor = ltmp::internal::thermal_conductivity/(ltmp::internal::micro_cell_size*ltmp::internal::micro_cell_size*1.e-20);

         // Determine change in Te and Tp
         for(int cell=0; cell<ltmp::internal::attenuation_array.size(); ++cell){

            const double Te = root_temperature_array[2*cell+0]*root_temperature_array[2*cell+0];
            const double Tp = root_temperature_array[2*cell+1]*root_temperature_array[2*cell+1];

            // calculate heat transfer from neighbouring cells
            double dTdiff = 0.0;
            for(int id=ltmp::internal::cell_neighbour_start_index[cell]; id<ltmp::internal::cell_neighbour_end_index[cell]; ++id){
               const int ncell = ltmp::internal::cell_neighbour_list[id]; // neighbour cell id
               double nTe = root_temperature_array[2*ncell+0]*root_temperature_array[2*ncell+0];
               dTdiff += nTe - Te;
            }

            delta_temperature_array[2*cell+0] = (G*(Tp-Te) + pump*attenuation_array[cell] + dTdiff*dTdiff_prefactor)*dt/(Ce*Te);
            delta_temperature_array[2*cell+1] = (G*(Te-Tp)                             )*dt/Cl;

         } // end of cell loop

         // Calculate new electron and lattice temperatures
         for(int cell=0; cell<ltmp::internal::attenuation_array.size(); ++cell){

            const double Te = root_temperature_array[2*cell+0]*root_temperature_array[2*cell+0] + delta_temperature_array[2*cell+0];
            const double Tp = root_temperature_array[2*cell+1]*root_temperature_array[2*cell+1] + delta_temperature_array[2*cell+1];

            root_temperature_array[2*cell+0] = sqrt(Te);
            root_temperature_array[2*cell+1] = sqrt(Tp);
         }

         // optionally output cell data
         if(ltmp::internal::output_microcell_data) ltmp::internal::write_cell_temperature_data();

         return;

      }

   } // end of namespace internal
} // end of namespace ltmp
