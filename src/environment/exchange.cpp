//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2016. All rights reserved.
//
//   Email: sj681@york.ac.uk
//
//------------------------------------------------------------------------------
//

// Vampire headers
#include "environment.hpp"
#include "internal.hpp"

// environment module headers
#include <stdlib.h>
#include <vector>
#include "errors.hpp"
#include "vio.hpp"
#include "sim.hpp"
#include "micromagnetic.hpp"
#include "../micromagnetic/internal.hpp"


namespace mm = micromagnetic::internal;


namespace environment
{

   namespace internal
   {

      //calcualtes the field from the micromagnetic cell onto the environment cell.
      std::vector < double > calculate_field_env(int env_cell, int mm_cell){

         double exchange_field[3]={0.0,0.0,0.0};
         //saves the x,y,z prefactors for the exchange constant
         double Ms_env = Ms;
         const double Acx = A*2/(Ms_env)*cell_size[2]*cell_size[1];
         const double Acy = A*2/(Ms_env)*cell_size[2]*cell_size[0];
         const double Acz = A*2/(Ms_env)*cell_size[0]*cell_size[1];

         //calculate |mj|
         const double mj = sqrt(mm::x_array[mm_cell]*mm::x_array[mm_cell] + mm::y_array[mm_cell]*mm::y_array[mm_cell] + mm::z_array[mm_cell]*mm::z_array[mm_cell]);
         //calcaultes the temperature dependant terms
         const double ECx = pow(mj,1.66)*Acx;
         const double ECy = pow(mj,1.66)*Acy;
         const double ECz = pow(mj,1.66)*Acz;
         //calcualtes the exchange field from i to j
         exchange_field[0] -= ECx*(mm::x_array[mm_cell] - x_array[env_cell]]);
         exchange_field[1] -= ECy*(mm::y_array[mm_cell] - y_array[env_cell]]);
         exchange_field[2] -= ECz*(mm::z_array[mm_cell] - z_array[env_cell]]);


         return exchange_field;
      }

      //calcualtes the field from the micromagnetic cell onto the environment cell.
      std::vector < double > calculate_field_mm(int env_cell, int mm_cell){

         double exchange_field[3]={0.0,0.0,0.0};
         //saves the x,y,z prefactors for the exchange constant
         const double Acx = A*2/(mm::Ms[cell])*mm::cell_size[2]*mm::cell_size[1];
         const double Acy = A*2/(mm::Ms[cell])*mm::cell_size[2]*mm::cell_size[0];
         const double Acz = A*2/(mm::Ms[cell])*mm::cell_size[0]*mm::cell_size[1];

         //calculate |mj|
         const double mj = sqrt(x_array[env_cell]*mm::x_array[env_cell] + mm::y_array[env_cell]*mm::y_array[env_cell] + mm::z_array[env_cell]*mm::z_array[env_cell]);
         //calcaultes the temperature dependant terms
         const double ECx = pow(mj,1.66)*Acx;
         const double ECy = pow(mj,1.66)*Acy;
         const double ECz = pow(mj,1.66)*Acz;
         //calcualtes the exchange field from i to j
         exchange_field[0] -= ECx*(x_array[env_cell] - mm::x_array[mm_cell]]);
         exchange_field[1] -= ECy*(y_array[env_cell] - mm::y_array[mm_cell]]);
         exchange_field[2] -= ECz*(z_array[env_cell] - mm::z_array[mm_cell]]);


         return exchange_field;
      }
   }
}
