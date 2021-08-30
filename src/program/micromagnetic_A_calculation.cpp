//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans and Andrea Meo 2014-2018. All rights reserved.
//
//-----------------------------------------------------------------------------
//

// Standard Libraries
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
#include "exchange.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"
#include "unitcell.hpp"
#include "../exchange/internal.hpp"
#include "../unitcell/internal.hpp"

namespace program{

//------------------------------------------------------------------------------
// Program to calculate a simple time series
//------------------------------------------------------------------------------
void mm_A_calculation(){

	// check calling of routine if error checking is activated
	if(err::check==true) std::cout << "program::A_calcultion has been called" << std::endl;

   double Jx,Jy,Jz;
   const double atomic_volume =  1e-30*unitcell::internal::unit_cell_size_x*unitcell::internal::unit_cell_size_x*unitcell::internal::unit_cell_size_x/cs::unit_cell.atom.size();
      double Ax = 0;
      double Ay = 0;
      double Az = 0;
   for(int atom = 0; atom < atoms::num_atoms; atom++){
   const int imaterial = atoms::type_array[atom];

      for(int nn = atoms::neighbour_list_start_index[atom];nn <= atoms::neighbour_list_end_index[atom]; nn++){
         const int natom = atoms::neighbour_list_array[nn];
         //const int jmaterial = atoms::type_array[natom];
         double dx = (atoms::x_coord_array[atom] - atoms::x_coord_array[natom]);
         double dy = (atoms::y_coord_array[atom] - atoms::y_coord_array[natom]);
         double dz = (atoms::z_coord_array[atom] - atoms::z_coord_array[natom]);
          if (dx <  -10)
             dx = cs::system_dimensions[0]  + dx;
         if (dy <  -10)
             dy = cs::system_dimensions[1]  + dy;
         if (dz <  -10)
             dz = cs::system_dimensions[2]  + dz;
         if (dx >  10)
             dx = cs::system_dimensions[0]  - dx;
         if (dy >  10)
             dy = cs::system_dimensions[1]  - dy;
         if (dz >  10)
             dz = cs::system_dimensions[2]  - dz;
         switch(exchange::internal::exchange_type){
            case exchange::isotropic:
               Jx = atoms::i_exchange_list[nn].Jij*mp::material[imaterial].mu_s_SI;
               Jy = atoms::i_exchange_list[nn].Jij*mp::material[imaterial].mu_s_SI;
               Jz = atoms::i_exchange_list[nn].Jij*mp::material[imaterial].mu_s_SI;

            break;
            case exchange::vectorial:
                  Jx = atoms::v_exchange_list[nn].Jij[0]*mp::material[imaterial].mu_s_SI;
                  Jy = atoms::v_exchange_list[nn].Jij[1]*mp::material[imaterial].mu_s_SI;
                  Jz = atoms::v_exchange_list[nn].Jij[2]*mp::material[imaterial].mu_s_SI;
            break;
            case exchange::tensorial:

                  Jx = atoms::t_exchange_list[nn].Jij[0][0]*mp::material[imaterial].mu_s_SI;
                  Jy = atoms::t_exchange_list[nn].Jij[1][1]*mp::material[imaterial].mu_s_SI;
                  Jz = atoms::t_exchange_list[nn].Jij[2][2]*mp::material[imaterial].mu_s_SI;
             break;
         }

         Ax = Ax + Jx*dx*dx*1e-20;
         Ay = Ay + Jy*dy*dy*1e-20;
         Az = Az + Jz*dz*dz*1e-20;
      }

   }

   std::cout <<  "Exchange constants\t Ax:\t"  << Ax/(atoms::num_atoms*atomic_volume) << "\tAy:\t" << Ay/(atoms::num_atoms*atomic_volume) << "\tAz:\t" << Az/(atoms::num_atoms*atomic_volume) << std::endl;
}

}
