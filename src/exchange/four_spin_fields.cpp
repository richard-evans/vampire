//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Mara Strungaru 2023. All rights reserved.
//
//   Email: mara.strungaru@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "atoms.hpp" // for exchange list type defs
#include "sim.hpp"
#include "exchange.hpp"

// exchange module headers
#include "internal.hpp"

double dot_product(double six, double siy, double siz, double sjx, double sjy, double sjz){
   double dp = six*sjx + siy*sjy + siz*sjz;
   return dp;
}

namespace exchange{

namespace internal{

void four_spin_exchange_fields(const int start_index, // first atom for exchange interactions to be calculated
                               const int end_index,                                 // field vectors for atoms
                               std::vector<double>& field_array_x,
                               std::vector<double>& field_array_y,
                               std::vector<double>& field_array_z){ // last +1 atom to be calculated){

   std::vector < int > numbers(atoms::num_atoms,0);

   // loop over all neighbours
	for(int nn = 0; nn < four_spin_neighbour_list_array_l.size(); ++nn){

      // get neighbouring atom number
   	const int atom = four_spin_neighbour_list_array_i[nn];
   	const int natomj = four_spin_neighbour_list_array_j[nn];
   	const int natomk = four_spin_neighbour_list_array_k[nn];
      const int natoml = four_spin_neighbour_list_array_l[nn];
      //const int jmaterial = atoms::type_array[natomj];
      const double Jij = four_spin_exchange_list[nn];

      const double sjx = atoms::x_spin_array[natomj];
      const double sjy = atoms::y_spin_array[natomj];
      const double sjz = atoms::z_spin_array[natomj];

      const double skx = atoms::x_spin_array[natomk];
      const double sky = atoms::y_spin_array[natomk];
      const double skz = atoms::z_spin_array[natomk];

      const double slx = atoms::x_spin_array[natoml];
      const double sly = atoms::y_spin_array[natoml];
      const double slz = atoms::z_spin_array[natoml];

      const double sk_dot_sl = dot_product(skx,sky,skz,slx,sly,slz);
      const double sj_dot_sk = dot_product(skx,sky,skz,sjx,sjy,sjz);
      const double sj_dot_sl = dot_product(sjx,sjy,sjz,slx,sly,slz);

      double athird=1.0/3.0;

      field_array_x[atom] = field_array_x[atom] + (Jij*athird)*(sjx*sk_dot_sl + skx*sj_dot_sl + slx*sj_dot_sk);
      field_array_y[atom] = field_array_y[atom] + (Jij*athird)*(sjy*sk_dot_sl + sky*sj_dot_sl + sly*sj_dot_sk);
      field_array_z[atom] = field_array_z[atom] + (Jij*athird)*(sjz*sk_dot_sl + skz*sj_dot_sl + slz*sj_dot_sk);

   }

   return;

}

} // end of internal namespace
} // end of exchange namespace
