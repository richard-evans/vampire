//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo, Sarah Jenkins and Richard F L Evans 2020. All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <cstdlib>
#include <iostream>

// Vampire headers
#include "dipole.hpp"
#include "vio.hpp"

// dipole module headers
#include "internal.hpp"

namespace dipole{

   namespace internal{

      //------------------------------------------------------------------------
      // Function to calculate Intra component of dipole tensors.
      //
      // The tensors between local cells with the cutoff range are calculated
      // explictly from the atomistic coordinates. Longer range tensors assume
      // the dipole-dipole form.
      //------------------------------------------------------------------------
      void compute_intra_tensor(const int celli,
                                const int cellj,
                                const int lc,
                                const std::vector<int>& global_atoms_in_cell_count,             // number of atoms in each cell (all CPUs)
                                const std::vector<int>& list_of_cells_with_atoms,               // list of cells to access atoms
                                const std::vector< std::vector<double> >& atoms_in_cells_array  // output array of positions and moments of atoms in cells
                               ){


         // Here we start
         // initialise temp vectors
         double tmp_rij_intra_xx = 0.0;
         double tmp_rij_intra_xy = 0.0;
         double tmp_rij_intra_xz = 0.0;

         double tmp_rij_intra_yy = 0.0;
         double tmp_rij_intra_yz = 0.0;
         double tmp_rij_intra_zz = 0.0;

         const int num_atoms = global_atoms_in_cell_count[celli];

         // search for cells i and j in local atom-cells list
         int cell_with_atoms_index_i = -1;
         for(size_t idx = 0; idx < atoms_in_cells_array.size(); idx++){
            const int cell = list_of_cells_with_atoms[idx];
            if( cell == celli ) cell_with_atoms_index_i = idx;
         }

         // check that proper cell is found
         if( cell_with_atoms_index_i == -1 ){
            std::cerr << "Programmer error! cell " << celli << " is not found in list of local cells with atomic positions!" << std::endl;
         }
         const int ci = cell_with_atoms_index_i;

         // loop over all atoms in cell i
         for(int pi = 0; pi < num_atoms; pi++){

            const double cix = atoms_in_cells_array[ci][4*pi+0];
            const double ciy = atoms_in_cells_array[ci][4*pi+1];
            const double ciz = atoms_in_cells_array[ci][4*pi+2];

            // loop over all atoms in cell for j < i (do half a full i-j loop)
            for( int qj = 0; qj < pi; qj++){

               const double rx = atoms_in_cells_array[ci][4*qj+0] - cix;
               const double ry = atoms_in_cells_array[ci][4*qj+1] - ciy;
               const double rz = atoms_in_cells_array[ci][4*qj+2] - ciz;

               const double rij = 1.0/sqrt(rx*rx+ry*ry+rz*rz); //Reciprocal of the distance
               const double rij3 = (rij*rij*rij); // Angstroms

               // for i -> j
               const double ex = rx*rij;
               const double ey = ry*rij;
               const double ez = rz*rij;

               // for j->i ex -> -ex, but squared in tensor so simple factor 2
               tmp_rij_intra_xx += 2.0*((3.0*ex*ex - 1.0)*rij3);
               tmp_rij_intra_xy += 2.0*((3.0*ex*ey      )*rij3);
               tmp_rij_intra_xz += 2.0*((3.0*ex*ez      )*rij3);

               tmp_rij_intra_yy += 2.0*((3.0*ey*ey - 1.0)*rij3);
               tmp_rij_intra_yz += 2.0*((3.0*ey*ez      )*rij3);
               tmp_rij_intra_zz += 2.0*((3.0*ez*ez - 1.0)*rij3);

            }
         }

         // normalisation factor accounting for i/j interactions (only symmetry of tensor is important)
         const double inorm = 1.0 / ( double(num_atoms) * double(num_atoms) );

         dipole::internal::rij_tensor_xx[lc][celli] =  (tmp_rij_intra_xx) * inorm;
         dipole::internal::rij_tensor_xy[lc][celli] =  (tmp_rij_intra_xy) * inorm;
         dipole::internal::rij_tensor_xz[lc][celli] =  (tmp_rij_intra_xz) * inorm;

         dipole::internal::rij_tensor_yy[lc][celli] =  (tmp_rij_intra_yy) * inorm;
         dipole::internal::rij_tensor_yz[lc][celli] =  (tmp_rij_intra_yz) * inorm;
         dipole::internal::rij_tensor_zz[lc][celli] =  (tmp_rij_intra_zz) * inorm;

         // Uncomment in case you want to check the tensor components
         // std::cout << "\n############# INTRA ###################\n";
         // std::cout << "lc = " << lc << "\ti = " << celli << std::endl;
         // std::cout << tmp_rij_intra_xx << "\t" << tmp_rij_intra_xy << "\t" << tmp_rij_intra_xz << "\n";
         // std::cout << tmp_rij_intra_xy << "\t" << tmp_rij_intra_yy << "\t" << tmp_rij_intra_yz << "\n";
         // std::cout << tmp_rij_intra_xz << "\t" << tmp_rij_intra_yz << "\t" << tmp_rij_intra_zz << "\n";
         // std::cout << "\n################################\n";
         // std::cout << std::endl;

      }  // End of function calculating Intra component of dipole tensor

   } // End of namespace internal
} // End of namespace dipole
