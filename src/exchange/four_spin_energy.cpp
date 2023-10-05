//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2018. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "atoms.hpp" // for exchange list type defs
#include "exchange.hpp"

// exchange module headers
#include "internal.hpp"

double dot_product2(double six, double siy, double siz, double sjx, double sjy, double sjz){
   double dp = six*sjx + siy*sjy + siz*sjz;
   return dp;
} 

namespace exchange{
   
   //-----------------------------------------------------------------------------------------
   // Function to calculate Four Spin exchange fields for spins between start and end index
   //-----------------------------------------------------------------------------------------
   // Four SPin exchange given by
   //
   //    E_4s = (Si.Sj)(Sk.Sl) + (Si.Sk)(Sj.Sl) + (Si.Sl)(Sj.Sk)
   //
   // Field given by:
   //
   //    H_4s^x = -dE/dSix =1/3 Dq (Sxj (Sk . Sl) + Sxk (Sl . Sj) + Sxl (Sk . Sj))
   //    H_4s^y = -dE/dSiy =1/3 Dq (Syj (Sk . Sl) + Syk (Sl . Sj) + Syl (Sk . Sj))
   //    H_4s^z = -dE/dSiz =1/3 Dq (Szj (Sk . Sl) + Szk (Sl . Sj) + Szl (Sk . Sj))
   //-----------------------------------------------------------------------------------------


   double single_spin_four_spin_energy(const int atom, const double sx, const double sy, const double sz){
         //std::cout << "here" <<std::endl;
// energy
   	double energy=0.0;

            const double six = atoms::x_spin_array[atom];
            const double siy = atoms::y_spin_array[atom];
            const double siz = atoms::z_spin_array[atom];

   	// Loop over neighbouring spins to calculate exchange
   	for(int nn = internal::four_spin_neighbour_list_start_index[atom]; nn <= internal::four_spin_neighbour_list_end_index[atom]; ++nn){

	         const int natomj = internal::four_spin_neighbour_list_array_j[nn];
				const int natomk = internal::four_spin_neighbour_list_array_k[nn];
				const int natoml = internal::four_spin_neighbour_list_array_l[nn];

            const double sjx = atoms::x_spin_array[natomj];
            const double sjy = atoms::y_spin_array[natomj];
            const double sjz = atoms::z_spin_array[natomj];

            const double skx = atoms::x_spin_array[natomk];
            const double sky = atoms::y_spin_array[natomk];
            const double skz = atoms::z_spin_array[natomk];

            const double slx = atoms::x_spin_array[natoml];
            const double sly = atoms::y_spin_array[natoml];
            const double slz = atoms::z_spin_array[natoml];

            const double si_dot_sj = dot_product2(six,siy,siz,sjx,sjy,sjz);
            const double si_dot_sk = dot_product2(six,siy,siz,skx,sky,skz);
            const double si_dot_sl = dot_product2(six,siy,siz,slx,sly,slz);
            const double sk_dot_sl = dot_product2(skx,sky,skz,slx,sly,slz);
            const double sj_dot_sk = dot_product2(skx,sky,skz,sjx,sjy,sjz);
            const double sj_dot_sl = dot_product2(sjx,sjy,sjz,slx,sly,slz);
            const double Jij = internal::four_spin_exchange_list[nn];
            energy = energy - 4.0*Jij/12.0*(si_dot_sj*sk_dot_sl + si_dot_sk*sj_dot_sl + si_dot_sl*sj_dot_sk);
       //     std::cout << energy << std::endl;
     }
      return 0;





}
}