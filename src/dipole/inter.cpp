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
      // Function to calculate inter component of dipole tensors.
      //
      // The tensors between local cells with the cutoff range are calculated
      // explictly from the atomistic coordinates. Longer range tensors assume
      // the dipole-dipole form.
      //------------------------------------------------------------------------
      void compute_inter_tensor(const int celli,                                                // global ID of cell i
                                const int cellj,                                                // global ID of cell i
                                const int lc,                                                   // local ID of cell
                                const double cutoff,                                            // cutoff range for dipole tensor construction (Angstroms)
                                const std::vector<int>& global_atoms_in_cell_count,             // number of atoms in each cell (all CPUs)
                                const std::vector<double>& cells_pos_and_mom_array,             // array of positions and cell moments
                                const std::vector<int>& list_of_cells_with_atoms,               // list of cells to access atoms
                                const std::vector< std::vector<double> >& atoms_in_cells_array  // output array of positions and moments of atoms in cells
                                ){

         // create temporary variables to store components of tensor
         double tmp_rij_inter_xx = 0.0;
         double tmp_rij_inter_xy = 0.0;
         double tmp_rij_inter_xz = 0.0;

         double tmp_rij_inter_yy = 0.0;
         double tmp_rij_inter_yz = 0.0;
         double tmp_rij_inter_zz = 0.0;

         // Calculate distance vectors between cells
         const double rx = cells_pos_and_mom_array[4*cellj+0] - cells_pos_and_mom_array[4*celli+0];
         const double ry = cells_pos_and_mom_array[4*cellj+1] - cells_pos_and_mom_array[4*celli+1];
         const double rz = cells_pos_and_mom_array[4*cellj+2] - cells_pos_and_mom_array[4*celli+2];

         // calculate square of distance for cutoff evaluation
         const double r2 = rx*rx + ry*ry + rz*rz;

         // If distance between macro-cells > cutoff nm => continuum approach (bare macro-cell method)
         if( r2 > cutoff*cutoff ){

            const double rij = 1.0/sqrt(r2); // Reciprocal of the distance

            //std::cout << dipole::atomistic_cutoff <<std::endl;
            // define unitarian distance vectors
	         const double ex = rx*rij;
	         const double ey = ry*rij;
	         const double ez = rz*rij;

	         const double rij3 = (rij*rij*rij); // Angstroms

            // calculate dipolar matrix for 6 entries because of symmetry
	         dipole::internal::rij_tensor_xx[lc][cellj] = ((3.0*ex*ex - 1.0)*rij3);
	         dipole::internal::rij_tensor_xy[lc][cellj] = ( 3.0*ex*ey      )*rij3 ;
	         dipole::internal::rij_tensor_xz[lc][cellj] = ( 3.0*ex*ez      )*rij3 ;

	         dipole::internal::rij_tensor_yy[lc][cellj] = ((3.0*ey*ey - 1.0)*rij3);
	         dipole::internal::rij_tensor_yz[lc][cellj] = ( 3.0*ey*ez      )*rij3 ;
	         dipole::internal::rij_tensor_zz[lc][cellj] = ((3.0*ez*ez - 1.0)*rij3);

         }

         //--------------------------------------------------------------------------
         // If distance between macro-cells < cutoff ==> apply inter-intra method
         //--------------------------------------------------------------------------
         else if( r2 <= cutoff*cutoff){

            const int num_i_atoms = global_atoms_in_cell_count[celli];
            const int num_j_atoms = global_atoms_in_cell_count[cellj];

            // search for cells i and j in local atom-cells list
            int cell_with_atoms_index_i = -1;
            int cell_with_atoms_index_j = -1;
            for(size_t idx = 0; idx < atoms_in_cells_array.size(); idx++){
               const int cell = list_of_cells_with_atoms[idx];
               if( cell == celli ) cell_with_atoms_index_i = idx;
               if( cell == cellj ) cell_with_atoms_index_j = idx;
            }

            // check that proper cell is found
            if( cell_with_atoms_index_i == -1 ){
               std::cerr << "Programmer error! cell " << celli << " is not found in list of local cells with atomic positions!" << std::endl;
            }
            if( cell_with_atoms_index_j == -1 ){
               std::cerr << "Programmer error! cell " << cellj << " is not found in list of local cells with atomic positions!" << std::endl;
            }
            const int ci = cell_with_atoms_index_i;
            const int cj = cell_with_atoms_index_j;

            // loop over all atoms in cell i
            for(int pi = 0; pi < num_i_atoms; pi++){

               const double cix = atoms_in_cells_array[ci][4*pi+0];
               const double ciy = atoms_in_cells_array[ci][4*pi+1];
               const double ciz = atoms_in_cells_array[ci][4*pi+2];

               // loop over all atoms in cell j
               for( int qj = 0; qj < num_j_atoms; qj++){

                  const double rx = atoms_in_cells_array[cj][4*qj+0] - cix;
                  const double ry = atoms_in_cells_array[cj][4*qj+1] - ciy;
                  const double rz = atoms_in_cells_array[cj][4*qj+2] - ciz;

                  const double dist = sqrt(rx*rx+ry*ry+rz*rz);
                  const double rij = 1.0/dist;  //Reciprocal of the distance

                  // if (dist <  dipole::atomistic_cutoff ){
                  // add to Jij somehow :P

                  const double ex = rx*rij;
                  const double ey = ry*rij;
                  const double ez = rz*rij;

                  const double rij3 = (rij*rij*rij); // Angstroms

                  tmp_rij_inter_xx += ((3.0*ex*ex - 1.0)*rij3);
                  tmp_rij_inter_xy += ((3.0*ex*ey      )*rij3);
                  tmp_rij_inter_xz += ((3.0*ex*ez      )*rij3);

                  tmp_rij_inter_yy += ((3.0*ey*ey - 1.0)*rij3);
                  tmp_rij_inter_yz += ((3.0*ey*ez      )*rij3);
                  tmp_rij_inter_zz += ((3.0*ez*ez - 1.0)*rij3);

               }
            }

            // normalisation factor accounting for i/j interactions (only symmetry of tensor is important)
            const double inorm = 1.0 / double( double(num_i_atoms) * double(num_j_atoms) );

            dipole::internal::rij_tensor_xx[lc][cellj] =  (tmp_rij_inter_xx) * inorm;
            dipole::internal::rij_tensor_xy[lc][cellj] =  (tmp_rij_inter_xy) * inorm;
            dipole::internal::rij_tensor_xz[lc][cellj] =  (tmp_rij_inter_xz) * inorm;

            dipole::internal::rij_tensor_yy[lc][cellj] =  (tmp_rij_inter_yy) * inorm;
            dipole::internal::rij_tensor_yz[lc][cellj] =  (tmp_rij_inter_yz) * inorm;
            dipole::internal::rij_tensor_zz[lc][cellj] =  (tmp_rij_inter_zz) * inorm;

            //if (i == 0) std::cout << "atom" <<  '\t' << i <<'\t' << j << "\t" << dipole::internal::rij_tensor_xx[lc][j] << "\t" << dipole::internal::rij_tensor_xy[lc][j] << '\t' <<dipole::internal::rij_tensor_xz[lc][j] << std::endl;
            // Uncomment in case you want to print the tensor components
            // std::cout << "\n############# INTER ###################\n";
            // std::cout << "lc = " << lc << "\tj = " << cellj << "\tNat_i\t" << num_i_atoms << "\tNat_j\t" << num_j_atoms << std::endl;
            // std::cout << tmp_rij_inter_xx << "\t" << tmp_rij_inter_xy << "\t" << tmp_rij_inter_xz << "\n";
            // std::cout << tmp_rij_inter_xy << "\t" << tmp_rij_inter_yy << "\t" << tmp_rij_inter_yz << "\n";
            // std::cout << tmp_rij_inter_xz << "\t" << tmp_rij_inter_yz << "\t" << tmp_rij_inter_zz << "\n";
            // std::cout << "\n################################\n";
            // std::cout << std::endl;

         }  // End of Inter part calculated atomicstically

            return;

      }  // End of funtion calculating inter component of dipole tensor

   } // End of namespace internal
} // End of namespace dipole
