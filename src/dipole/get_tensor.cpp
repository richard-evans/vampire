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
   std::vector<double> get_tensor_1D_xx(){
      std::vector<double> tensor_unrolled_xx;
      for(int lc=0; lc<dipole::internal::cells_num_local_cells; lc++){
      	for(unsigned int j=0; j<dipole::internal::rij_tensor_xx[lc].size(); j++){
				tensor_unrolled_xx.push_back( dipole::internal::rij_tensor_xx[lc][j] );
      	}
      }
      return tensor_unrolled_xx;
   }

   std::vector<double> get_tensor_1D_xy(){
      std::vector<double> tensor_unrolled_xy;
      for(int lc=0; lc<dipole::internal::cells_num_local_cells; lc++){
      	for(unsigned int j=0; j<dipole::internal::rij_tensor_xy[lc].size(); j++){
				tensor_unrolled_xy.push_back( dipole::internal::rij_tensor_xy[lc][j] );
      	}
      }
      return tensor_unrolled_xy;
   }

   std::vector<double> get_tensor_1D_xz(){
      std::vector<double> tensor_unrolled_xz;
      for(int lc=0; lc<dipole::internal::cells_num_local_cells; lc++){
      	for(unsigned int j=0; j<dipole::internal::rij_tensor_xz[lc].size(); j++){
				tensor_unrolled_xz.push_back( dipole::internal::rij_tensor_xz[lc][j] );
      	}
      }
      return tensor_unrolled_xz;
   }

   std::vector<double> get_tensor_1D_yy(){
      std::vector<double> tensor_unrolled_yy;
      for(int lc=0; lc<dipole::internal::cells_num_local_cells; lc++){
      	for(unsigned int j=0; j<dipole::internal::rij_tensor_yy[lc].size(); j++){
				tensor_unrolled_yy.push_back( dipole::internal::rij_tensor_yy[lc][j] );
      	}
      }
      return tensor_unrolled_yy;
   }

   std::vector<double> get_tensor_1D_yz(){
      std::vector<double> tensor_unrolled_yz;
      for(int lc=0; lc<dipole::internal::cells_num_local_cells; lc++){
      	for(unsigned int j=0; j<dipole::internal::rij_tensor_yz[lc].size(); j++){
				tensor_unrolled_yz.push_back( dipole::internal::rij_tensor_yz[lc][j] );
      	}
      }
      return tensor_unrolled_yz;
   }

   std::vector<double> get_tensor_1D_zz(){
      std::vector<double> tensor_unrolled_zz;
      for(int lc=0; lc<dipole::internal::cells_num_local_cells; lc++){
      	for(unsigned int j=0; j<dipole::internal::rij_tensor_zz[lc].size(); j++){
				tensor_unrolled_zz.push_back( dipole::internal::rij_tensor_zz[lc][j] );
      	}
      }
      return tensor_unrolled_zz;
   }

} // end of dipole namespace




