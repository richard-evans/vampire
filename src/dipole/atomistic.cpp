//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2019. All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <iostream>

// Vampire headers
#include "dipole.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "material.hpp"

// dipole module headers
#include "internal.hpp"

// alias interal dipole namespace for brevity
namespace dp = dipole::internal;

namespace dipole{

   namespace internal{

      //------------------------------------------------------------------------
      // Function to initialise atomistic resolution dipole solver
      //------------------------------------------------------------------------
      void initialize_atomistic_solver(int num_atoms,                      // number of atoms (only correct in serial)
                                       std::vector<double>& x_coord_array, // atomic corrdinates (angstroms)
                                       std::vector<double>& y_coord_array,
                                       std::vector<double>& z_coord_array,
                                       std::vector<double>& moments_array, // atomistic magnetic moments (bohr magnetons)
                                       std::vector<int>& mat_id_array){    // atom material ID

         // Calculate memory requirements and inform user
         const double mem = double(num_atoms) * double(vmpi::num_processors) * sizeof(double) * 7.0 / 1.0e6;
         zlog << zTs() << "Atomistic dipole field calculation has been enabled and requires " << mem << " MB of RAM" << std::endl;
         std::cout     << "Atomistic dipole field calculation has been enabled and requires " << mem << " MB of RAM" << std::endl;

      #ifdef MPICF

         // calculate number of local atoms excluding halo
         dp::num_local_atoms = vmpi::num_core_atoms + vmpi::num_bdry_atoms;

         // Zero all moments for non-magnetic materials by using copy of moments array
         std::vector<double> moments_array_copy = moments_array;
         for (int atom = 0; atom < dp::num_local_atoms; atom++){
            int mat = mat_id_array[atom]; // get atomic material ID
            if(mp::material[mat].non_magnetic != 0) moments_array_copy[atom] = 0.0; // set local moment to zero
         }

         // resize main displacement array
         dp::receive_displacements.resize(vmpi::num_processors,0);
         dp::receive_counts.resize(vmpi::num_processors,0);

         // Collate the number of atoms from each process on root
         MPI_Gather(&dp::num_local_atoms, 1, MPI_INT, &dp::receive_counts[0], 1, MPI_INT, 0, MPI_COMM_WORLD);

         // calculate the total number of atoms sent
         for (int proc = 1; proc < vmpi::num_processors; proc ++){
            dp::receive_displacements[proc] = dp::receive_displacements[proc-1] + dp::receive_counts[proc-1];
         }

         // Broadcast displacements and counts to all processors
         MPI_Bcast(&dp::receive_displacements[0], vmpi::num_processors, MPI_INT, 0, MPI_COMM_WORLD);
         MPI_Bcast(&dp::receive_counts[0],        vmpi::num_processors, MPI_INT, 0, MPI_COMM_WORLD);

         // calculate total number of atoms from displacements
         dp::total_num_atoms = dp::receive_displacements[vmpi::num_processors-1] + dp::receive_counts[vmpi::num_processors - 1];

         // broadcast total number of atoms to all processes
         MPI_Bcast(&dp::total_num_atoms, 1, MPI_INT, 0, MPI_COMM_WORLD);

         //std::cerr << vmpi::my_rank << "\t" << total_num_atoms << "\t" << num_local_atoms << std::endl;

         // resize coordinate arrays to include all atoms
         dp::cx.resize(total_num_atoms, 0.0);
         dp::cy.resize(total_num_atoms, 0.0);
         dp::cz.resize(total_num_atoms, 0.0);
         dp::sm.resize(total_num_atoms, 0.0);

         // Gather atomic positions on all processors
         MPI_Allgatherv(&x_coord_array[0],      num_local_atoms, MPI_DOUBLE, &dp::cx[0], &dp::receive_counts[0], &dp::receive_displacements[0], MPI_DOUBLE, MPI_COMM_WORLD);
         MPI_Allgatherv(&y_coord_array[0],      num_local_atoms, MPI_DOUBLE, &dp::cy[0], &dp::receive_counts[0], &dp::receive_displacements[0], MPI_DOUBLE, MPI_COMM_WORLD);
         MPI_Allgatherv(&z_coord_array[0],      num_local_atoms, MPI_DOUBLE, &dp::cz[0], &dp::receive_counts[0], &dp::receive_displacements[0], MPI_DOUBLE, MPI_COMM_WORLD);
         MPI_Allgatherv(&moments_array_copy[0], num_local_atoms, MPI_DOUBLE, &dp::sm[0], &dp::receive_counts[0], &dp::receive_displacements[0], MPI_DOUBLE, MPI_COMM_WORLD);

         // Resize arrays to hold all spin and moment positions
         dp::sx.resize(total_num_atoms, 0.0);
         dp::sy.resize(total_num_atoms, 0.0);
         dp::sz.resize(total_num_atoms, 1.0);

      #else

         // calculate number of local atoms
         dp::num_local_atoms = num_atoms;
         dp::total_num_atoms = num_local_atoms;

         // Zero all moments for non-magnetic materials by using copy of moments array
         std::vector<double> moments_array_copy = moments_array;
         for (int atom = 0; atom < dp::num_local_atoms; atom++){
            int mat = mat_id_array[atom]; // get atomic material ID
            if(mp::material[mat].non_magnetic != 0) moments_array_copy[atom] = 0.0; // set local moment to zero
         }

         // resize coordinate arrays to include all atoms
         dp::cx.resize(total_num_atoms, 0.0);
         dp::cy.resize(total_num_atoms, 0.0);
         dp::cz.resize(total_num_atoms, 0.0);
         dp::sm.resize(total_num_atoms, 0.0);

         // copy atomic coordinates to local array
         for (int atom = 0; atom < dp::total_num_atoms; atom++){
            cx[atom] = x_coord_array[atom];
            cy[atom] = y_coord_array[atom];
            cz[atom] = z_coord_array[atom];
            sm[atom] = moments_array_copy[atom];
         }

         // Resize arrays to hold all spin and moment positions
         dp::sx.resize(total_num_atoms, 0.0);
         dp::sy.resize(total_num_atoms, 0.0);
         dp::sz.resize(total_num_atoms, 1.0);

      #endif

      // If enabled, output calculated atomistic dipole field coordinates and moments (passing local values)
      if(dp::output_atomistic_dipole_field) output_atomistic_coordinates(num_atoms, x_coord_array, y_coord_array, z_coord_array, moments_array_copy);

      // Inform user of successful initialisation
      zlog << zTs() << "Atomistic dipole field calculation initialised" << std::endl;
      std::cout     << "Atomistic dipole field calculation initialised" << std::endl;

      return;

   }

   //------------------------------------------------------------------------
   // Function to calculate atomistic resolution dipole field
   //
   //
   //------------------------------------------------------------------------
   void calculate_atomistic_dipole_field(std::vector<double>& x_spin_array, // atomic spin directions
                                         std::vector<double>& y_spin_array,
                                         std::vector<double>& z_spin_array){

       const double prefactor = 0.9274009994; // mu_o_4pi * muB / Angstrom^3 = 1.0e-7 * 9.274009994e-24 / 1.0e-30 = 0.9274009994

       // cast number of local atoms to a local constant
       const int num_atoms_on_my_processor = dp::num_local_atoms;

      #ifdef MPICF

         // collate and broadcast new spin positions to all processors
         MPI_Allgatherv(&x_spin_array[0], num_local_atoms, MPI_DOUBLE, &dp::sx[0], &dp::receive_counts[0], &dp::receive_displacements[0], MPI_DOUBLE, MPI_COMM_WORLD);
         MPI_Allgatherv(&y_spin_array[0], num_local_atoms, MPI_DOUBLE, &dp::sy[0], &dp::receive_counts[0], &dp::receive_displacements[0], MPI_DOUBLE, MPI_COMM_WORLD);
         MPI_Allgatherv(&z_spin_array[0], num_local_atoms, MPI_DOUBLE, &dp::sz[0], &dp::receive_counts[0], &dp::receive_displacements[0], MPI_DOUBLE, MPI_COMM_WORLD);

      #else

          // copy spin directions to local array
          for (int atom = 0; atom < dp::total_num_atoms; atom++){
             dp::sx[atom] = x_spin_array[atom];
             dp::sy[atom] = y_spin_array[atom];
             dp::sz[atom] = z_spin_array[atom];
          }

      #endif

      //------------------------------------------------------------------------
      // Loop over all local atoms
      //------------------------------------------------------------------------
      for (int atom_i = 0; atom_i < num_atoms_on_my_processor; atom_i ++){

         // get id of atom i in total list
         #ifdef MPICF
          int atom_i_in_total_list_j = receive_displacements[vmpi::my_rank] + atom_i;
         #else
          int atom_i_in_total_list_j = atom_i;
         #endif

         // get coordinates of atom i
         const double xi = dp::cx[atom_i_in_total_list_j];
         const double yi = dp::cy[atom_i_in_total_list_j];
         const double zi = dp::cz[atom_i_in_total_list_j];

         // temporary variables to calculate local field
         double bx = 0.0;
         double by = 0.0;
         double bz = 0.0;

         //------------------------------------------------------------------------
         // Loop over all other atoms j < i
         //------------------------------------------------------------------------
         for (int atom_j = 0; atom_j < atom_i_in_total_list_j; atom_j ++){

            // get spin moment of atom j
            const double mj = dp::sm[atom_j];

            // get coordinates of atom j
            const double xj = dp::cx[atom_j];
            const double yj = dp::cy[atom_j];
            const double zj = dp::cz[atom_j];

            // calculate net spin moment of atom j
            const double mxj = dp::sx[atom_j] * mj;
            const double myj = dp::sy[atom_j] * mj;
            const double mzj = dp::sz[atom_j] * mj;

            // calculate position vector i -> j
            const double rx = xj - xi;
            const double ry = yj - yi;
            const double rz = zj - zi;

            // calculate distance between atoms
            const double rij = 1.0/sqrt(rx*rx+ry*ry+rz*rz); //Reciprocal of the distance

            // calculate unit vector from i -> j
            const double ex = rx*rij;
            const double ey = ry*rij;
            const double ez = rz*rij;

            // calculate cube of distance for normalisation
            const double rij3 = ( rij * rij * rij); // Angstroms

            // calculate r . m
            const double rdotm = ex*mxj + ey*myj + ez*mzj;

            bx += (3.0*ex*rdotm - mxj) * rij3;
            by += (3.0*ey*rdotm - myj) * rij3;
            bz += (3.0*ez*rdotm - mzj) * rij3;

         }

         //------------------------------------------------------------------------
         // Loop over all other atoms j > i
         //------------------------------------------------------------------------
         for (int atom_j = atom_i_in_total_list_j + 1; atom_j < dp::total_num_atoms; atom_j ++){

            // get spin moment of atom j
            const double mj = dp::sm[atom_j];

            // get coordinates of atom j
            const double xj = dp::cx[atom_j];
            const double yj = dp::cy[atom_j];
            const double zj = dp::cz[atom_j];

            // calculate net spin moment of atom j
            const double mxj = dp::sx[atom_j] * mj;
            const double myj = dp::sy[atom_j] * mj;
            const double mzj = dp::sz[atom_j] * mj;

            // calculate position vector i -> j
            const double rx = xj - xi;
            const double ry = yj - yi;
            const double rz = zj - zi;

            // calculate distance between atoms
            const double rij = 1.0/sqrt(rx*rx+ry*ry+rz*rz); //Reciprocal of the distance

            // calculate unit vector from i -> j
            const double ex = rx*rij;
            const double ey = ry*rij;
            const double ez = rz*rij;

            // calculate cube of distance for normalisation
            const double rij3 = ( rij * rij * rij); // Angstroms

            // calculate r . m
            const double rdotm = ex*mxj + ey*myj + ez*mzj;

            bx += (3.0*ex*rdotm - mxj) * rij3;
            by += (3.0*ey*rdotm - myj) * rij3;
            bz += (3.0*ez*rdotm - mzj) * rij3;

         }

         // save total dipole field to atomic field array
         dipole::atom_dipolar_field_array_x[atom_i] = prefactor * bx;
         dipole::atom_dipolar_field_array_y[atom_i] = prefactor * by;
         dipole::atom_dipolar_field_array_z[atom_i] = prefactor * bz;

      }

      if(dp::output_atomistic_dipole_field) output_atomistic_dipole_fields();

      return;

   }

} // end of namespace internal
} // end of namespace dipole
