//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) Andrea Meo Meo 2021. All rights reserved.
//
//-----------------------------------------------------------------------------


// C++ standard library headers
#include <iostream>

// Vampire headers
#include "dipole.hpp"

// exchange module headers
#include "internal.hpp"

namespace dipole{

   //------------------------------------------------------------------------------
   // Functions to return unrolled dipole tensor
   //------------------------------------------------------------------------------
   std::vector<double> unroll_tensor(const int element, double dummy){
      std::vector<double> out;
      // Select which component of tensor to unrol since it belongs to dipole::internal
      std::vector<std::vector<double> > in;
      switch(element){
         case 1:
            in = dipole::internal::rij_tensor_xx;
            break;
         case 2:
            in = dipole::internal::rij_tensor_xy;
            break;
         case 3:
            in = dipole::internal::rij_tensor_xz;
            break;
         case 4:
            in = dipole::internal::rij_tensor_yy;
            break;
         case 5:
            in = dipole::internal::rij_tensor_yz;
            break;
         case 6:
            in = dipole::internal::rij_tensor_zz;
            break;
      }
      for(int64_t lc=0; lc<dipole::internal::cells_num_local_cells; lc++){
      	for(size_t j=0; j<in[lc].size(); j++){ out.push_back( in[lc][j] ); }
      }
      return out;
   }

   std::vector<float> unroll_tensor(const int element, float dummy){
      std::vector<float> out;
      // Select which component of tensor to unrol since it belongs to dipole::internal
      std::vector<std::vector<double> > in;
      switch(element){
         case 1:
            in = dipole::internal::rij_tensor_xx;
            break;
         case 2:
            in = dipole::internal::rij_tensor_xy;
            break;
         case 3:
            in = dipole::internal::rij_tensor_xz;
            break;
         case 4:
            in = dipole::internal::rij_tensor_yy;
            break;
         case 5:
            in = dipole::internal::rij_tensor_yz;
            break;
         case 6:
            in = dipole::internal::rij_tensor_zz;
            break;
      }
      for(int64_t lc=0; lc<dipole::internal::cells_num_local_cells; lc++){
      	for(size_t j=0; j<in[lc].size(); j++){ out.push_back( in[lc][j] ); }
      }
      return out;
   }

}
