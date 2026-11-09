//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins 2018. All rights reserved.
//
//   Email: sarah.jenkins@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
// C++ standard library headers
#include <string>
#include <cmath>
#include <cstdlib>
#include <iostream>
// C library headers
#include <fenv.h>
#include <signal.h>
#include <math.h>
// Vampire headers
#include "cells.hpp" // needed for cells::cell_id_array but to be removed
#include "dipole.hpp"
#include "errors.hpp"
#include "vio.hpp"
#include "vutil.hpp"
#include "atoms.hpp"
// Vampire headers
#include "hierarchical.hpp"


// hierarchical module headers
#include "internal.hpp"

namespace ha = hierarchical::internal;

namespace hierarchical{

   //------------------------------------------------------------------------
   // Function to calculate Intra component of dipole tensors.
   //
   // The tensors between local cells with the cutoff range are calculated
   // explictly from the atomistic coordinates. Longer range tensors assume
   // the dipole-dipole form.
   //------------------------------------------------------------------------
  void ha::calc_intra(const int celli,
                      const int cellj,
                      const int interaction_no,
                      const std::vector<int>& global_atoms_in_cell_count,             // number of atoms in each cell (all CPUs)
                      const region_atoms_t& atoms_table                               // atoms of the cells near this process's cells
                     ){

      // initialise temp vectors
      double tmp_rij_intra_xx = 0.0;
      double tmp_rij_intra_xy = 0.0;
      double tmp_rij_intra_xz = 0.0;

      double tmp_rij_intra_yy = 0.0;
      double tmp_rij_intra_yz = 0.0;
      double tmp_rij_intra_zz = 0.0;

      const int num_atoms = global_atoms_in_cell_count[celli];

      // find cell i in the atom table of the region
      const int ki = atoms_table.find(celli);
      if( ki == -1 || atoms_table.cell_starts[ki+1] - atoms_table.cell_starts[ki] != num_atoms ){
         err::fatal("Programmer error: cell " + std::to_string(celli) + " is not complete in the hierarchical atom table");
      }
      const double* ci = &atoms_table.coords[3*atoms_table.cell_starts[ki]]; // coordinates of the atoms in cell i

      // loop over all atoms in cell i
      for(int pi = 0; pi < num_atoms; pi++){

         const double cix = ci[3*pi+0];
         const double ciy = ci[3*pi+1];
         const double ciz = ci[3*pi+2];

         // loop over all atoms in cell for j < i (do half a full i-j loop)
         for( int qj = 0; qj < pi; qj++){

            const double rx = ci[3*qj+0] - cix;
            const double ry = ci[3*qj+1] - ciy;
            const double rz = ci[3*qj+2] - ciz;

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

      ha::rij_tensor_xx[interaction_no] =  (tmp_rij_intra_xx) * inorm;
      ha::rij_tensor_xy[interaction_no] =  (tmp_rij_intra_xy) * inorm;
      ha::rij_tensor_xz[interaction_no] =  (tmp_rij_intra_xz) * inorm;

      ha::rij_tensor_yy[interaction_no] =  (tmp_rij_intra_yy) * inorm;
      ha::rij_tensor_yz[interaction_no] =  (tmp_rij_intra_yz) * inorm;
      ha::rij_tensor_zz[interaction_no] =  (tmp_rij_intra_zz) * inorm;

      //Uncomment in case you want to check the tensor components
      // std::cout << "\n############# INTRA ###################\n";
      // std::cout << "interaction = " << interaction_no << "\ti = " << celli << std::endl;
      // std::cout << tmp_rij_intra_xx << "\t" << tmp_rij_intra_xy << "\t" << tmp_rij_intra_xz << "\n";
      // std::cout << tmp_rij_intra_xy << "\t" << tmp_rij_intra_yy << "\t" << tmp_rij_intra_yz << "\n";
      // std::cout << tmp_rij_intra_xz << "\t" << tmp_rij_intra_yz << "\t" << tmp_rij_intra_zz << "\n";
      // std::cout << "\n################################\n";
      // std::cout << std::endl;

      return;

   }  // End of function calculating Intra component of dipole tensor

} // end of hierarchical namespace
