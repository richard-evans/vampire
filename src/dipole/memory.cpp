//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo and Richard F L Evans 2016. All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "dipole.hpp"

// dipole module headers
#include "internal.hpp"

namespace dipole{

   namespace internal{

      //-----------------------------------------------------------------
      // Function to allocate internal memory for dipole tensors
      //-----------------------------------------------------------------
      void allocate_memory(const int cells_num_local_cells, const int cells_num_cells){

         // reserve memory for inter cell arrays
         dipole::internal::rij_tensor_xx.reserve(cells_num_local_cells);
         dipole::internal::rij_tensor_xy.reserve(cells_num_local_cells);
         dipole::internal::rij_tensor_xz.reserve(cells_num_local_cells);
         dipole::internal::rij_tensor_yy.reserve(cells_num_local_cells);
         dipole::internal::rij_tensor_yz.reserve(cells_num_local_cells);
         dipole::internal::rij_tensor_zz.reserve(cells_num_local_cells);


         // allocate arrays to store data [nloccell x ncells]
         for(int lc=0; lc<cells_num_local_cells; lc++){

            dipole::internal::rij_tensor_xx.push_back(std::vector<double>());
            dipole::internal::rij_tensor_xx[lc].resize(cells_num_cells,0.0);

            dipole::internal::rij_tensor_xy.push_back(std::vector<double>());
            dipole::internal::rij_tensor_xy[lc].resize(cells_num_cells,0.0);

            dipole::internal::rij_tensor_xz.push_back(std::vector<double>());
            dipole::internal::rij_tensor_xz[lc].resize(cells_num_cells,0.0);

            dipole::internal::rij_tensor_yy.push_back(std::vector<double>());
            dipole::internal::rij_tensor_yy[lc].resize(cells_num_cells,0.0);

            dipole::internal::rij_tensor_yz.push_back(std::vector<double>());
            dipole::internal::rij_tensor_yz[lc].resize(cells_num_cells,0.0);

            dipole::internal::rij_tensor_zz.push_back(std::vector<double>());
            dipole::internal::rij_tensor_zz[lc].resize(cells_num_cells,0.0);
         }

         // resize B-field cells array
         dipole::cells_field_array_x.resize(cells_num_cells,0.0);
         dipole::cells_field_array_y.resize(cells_num_cells,0.0);
         dipole::cells_field_array_z.resize(cells_num_cells,0.0);

         // resize mu_0*Hd-field cells array
         dipole::cells_mu0Hd_field_array_x.resize(cells_num_cells,0.0);
         dipole::cells_mu0Hd_field_array_y.resize(cells_num_cells,0.0);
         dipole::cells_mu0Hd_field_array_z.resize(cells_num_cells,0.0);

         return;

      }

   } // end of namespace internal

} // end of namespace dipole
