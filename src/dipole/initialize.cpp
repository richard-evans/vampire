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
#include <cmath>
#include <cstdlib>
#include <iostream>

// Vampire headers
#include "cells.hpp"
#include "dipole.hpp"
#include "vio.hpp"
#include "vutil.hpp"


// dipole module headers
#include "internal.hpp"

namespace dipole{

   //----------------------------------------------------------------------------
   // Function to initialize dipole module
   //----------------------------------------------------------------------------
   void initialize(const int cells_num_atoms_in_unit_cell,
                  int cells_num_cells, /// number of macrocells
                  int cells_num_local_cells, /// number of local macrocells
                  const double cells_macro_cell_size,
                  std::vector <int>& cells_local_cell_array,
                  std::vector <int>& cells_num_atoms_in_cell, /// number of atoms in each cell
                  std::vector <int>& cells_num_atoms_in_cell_global, /// number of atoms in each cell
                  std::vector < std::vector <int> >& cells_index_atoms_array,
                  std::vector<double>& cells_volume_array,
                  std::vector<double>& cells_pos_and_mom_array, // array to store positions and moment of cells
                  std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_x,
                  std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_y,
                  std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_z,
                  std::vector<int>& atom_type_array,
                  std::vector<int>& atom_cell_id_array,
                  std::vector<double>& atom_coords_x, //atomic coordinates
                  std::vector<double>& atom_coords_y,
                  std::vector<double>& atom_coords_z,
                  std::vector<double>& x_spin_array, // atomic spin directions
                  std::vector<double>& y_spin_array,
                  std::vector<double>& z_spin_array,
                  std::vector<double>& atom_moments, // atomic magnetic moments
                  int num_atoms
				){

	   //-------------------------------------------------------------------------------------
      // Check for dipole calculation enabled, if not do nothing
      //-------------------------------------------------------------------------------------
      if(!dipole::activated) return;

		// check for prior initialisation
		if(dipole::internal::initialised){
      	zlog << zTs() << "Warning:  Dipole field calculation already initialised. Continuing." << std::endl;
      	return;
		}

      // output informative message
      std::cout << "Initialising dipole field calculation" << std::endl;
		zlog << zTs() << "Initialising dipole field calculation" << std::endl;

      // allocate memory for rij matrix
      dipole::internal::allocate_memory(cells_num_local_cells, cells_num_cells);

      //-------------------------------------------------------------------------------------
      // Set const for functions
      //-------------------------------------------------------------------------------------

      dipole::internal::num_atoms                  = num_atoms;
      dipole::internal::atom_type_array            = atom_type_array;
      dipole::internal::atom_cell_id_array         = atom_cell_id_array;

      dipole::internal::cells_num_cells            = cells_num_cells;
      dipole::internal::cells_num_local_cells      = cells_num_local_cells;
      dipole::internal::cells_local_cell_array     = cells_local_cell_array;
      dipole::internal::cells_num_atoms_in_cell    = cells_num_atoms_in_cell;
      dipole::internal::cells_volume_array         = cells_volume_array;

      dipole::internal::cells_pos_and_mom_array    = cells_pos_and_mom_array;

		//-------------------------------------------------------------------------------------
		// Starting calculation of dipolar field
		//-------------------------------------------------------------------------------------

      // Check memory requirements and print to screen
      zlog << zTs() << "Fast dipole field calculation has been enabled and requires " << double(dipole::internal::cells_num_cells)*double(dipole::internal::cells_num_local_cells*6)*8.0/1.0e6 << " MB of RAM" << std::endl;
      std::cout     << "Fast dipole field calculation has been enabled and requires " << double(dipole::internal::cells_num_cells)*double(dipole::internal::cells_num_local_cells*6)*8.0/1.0e6 << " MB of RAM" << std::endl;

      zlog << zTs() << "Total memory for dipole calculation (all CPUs): " << double(dipole::internal::cells_num_cells)*double(dipole::internal::cells_num_cells*6)*8.0/1.0e6 << " MB of RAM" << std::endl;
      std::cout << "Total memory for dipole calculation (all CPUs): " << double(dipole::internal::cells_num_cells)*double(dipole::internal::cells_num_cells*6)*8.0/1.0e6 << " MB of RAM" << std::endl;

      zlog << zTs() << "Number of local cells for dipole calculation = " << dipole::internal::cells_num_local_cells << std::endl;
      zlog << zTs() << "Number of total cells for dipole calculation = " << dipole::internal::cells_num_cells << std::endl;

      //----------------------------------------------------------
      // Calculation of dipolar tensor
      //----------------------------------------------------------
      switch (dipole::internal::solver){

         case dipole::internal::macrocell:
            dipole::internal::initialize_macrocell_solver(cells_num_atoms_in_unit_cell, dipole::internal::cells_num_cells, dipole::internal::cells_num_local_cells, cells_macro_cell_size, dipole::internal::cells_local_cell_array,
                                                       dipole::internal::cells_num_atoms_in_cell, cells_num_atoms_in_cell_global, cells_index_atoms_array, dipole::internal::cells_volume_array, dipole::internal::cells_pos_and_mom_array,
                                                       cells_atom_in_cell_coords_array_x, cells_atom_in_cell_coords_array_y, cells_atom_in_cell_coords_array_z,
                                                       dipole::internal::atom_type_array, dipole::internal::atom_cell_id_array, atom_coords_x, atom_coords_y, atom_coords_z, dipole::internal::num_atoms);
            break;

         case dipole::internal::tensor:
            dipole::internal::initialize_tensor_solver(cells_num_atoms_in_unit_cell, dipole::internal::cells_num_cells, dipole::internal::cells_num_local_cells, cells_macro_cell_size, dipole::internal::cells_local_cell_array,
                                                       dipole::internal::cells_num_atoms_in_cell, cells_num_atoms_in_cell_global, cells_index_atoms_array, dipole::internal::cells_volume_array, dipole::internal::cells_pos_and_mom_array,
                                                       cells_atom_in_cell_coords_array_x, cells_atom_in_cell_coords_array_y, cells_atom_in_cell_coords_array_z,
                                                       dipole::internal::atom_type_array, dipole::internal::atom_cell_id_array, atom_coords_x, atom_coords_y, atom_coords_z, dipole::internal::num_atoms);
            break;

         case dipole::internal::atomistic:
            dipole::internal::initialize_atomistic_solver(num_atoms, atom_coords_x, atom_coords_y, atom_coords_z, atom_moments, atom_type_array);
            break;

      }

      // Set initialised flag
      dipole::internal::initialised=true;

      //------------------------------------------------------------------------
      // Precalculate dipole field and time for performance
      //------------------------------------------------------------------------

      // instantiate timer
      vutil::vtimer_t timer;

      // start timer
      timer.start();

      // now calculate fields at zero time
      dipole::calculate_field(0, x_spin_array, y_spin_array, z_spin_array);

      // hold parallel calculation until all processors have completed the update
      vmpi::barrier();

      // stop timer
      timer.stop();

      zlog << zTs() << "Time required for dipole update: " << timer.elapsed_time() << " s." << std::endl;

      //--------------------------------------------------------------------------------------------------
      // Calculate gloabl demagnetizing factor from dipole tensors
      //--------------------------------------------------------------------------------------------------
      //
      //
      //--------------------------------------------------------------------------------------------------

      // Output informative message to log file
      zlog << zTs() << "Outputting dipole matrix " << std::endl;

      // Output Demag tensor only if first step of simulation since depending only on shape
      std::vector<double> N_tensor_array(6*dipole::internal::cells_num_cells,0.0);


      // Every cpus print to check dipolar matrix inter term
      for(int lc=0; lc<dipole::internal::cells_num_local_cells; lc++){

         // get id of cell
         int i = cells::cell_id_array[lc];

         // if local cell contains atoms
         if(dipole::internal::cells_num_atoms_in_cell[i]>0){

            // loop over all neighbours for cell
            for(unsigned int j=0; j<dipole::internal::rij_tensor_xx[lc].size(); j++){
               if(dipole::internal::cells_num_atoms_in_cell[j]>0){

                  // To obtain dipolar matrix free of units, multiply tensor by "factor"
                  //const double Vatomic = dipole::internal::cells_volume_array[j]/double(dipole::internal::cells_num_atoms_in_cell[j]);
                  const double factor = double(dipole::internal::cells_num_atoms_in_cell[j]) * double(dipole::internal::cells_num_atoms_in_cell[i]);

                  N_tensor_array[6*i+0] +=  factor*(dipole::internal::rij_tensor_xx[lc][j]);
                  N_tensor_array[6*i+1] +=  factor*(dipole::internal::rij_tensor_xy[lc][j]);
                  N_tensor_array[6*i+2] +=  factor*(dipole::internal::rij_tensor_xz[lc][j]);
                  N_tensor_array[6*i+3] +=  factor*(dipole::internal::rij_tensor_yy[lc][j]);
                  N_tensor_array[6*i+4] +=  factor*(dipole::internal::rij_tensor_yz[lc][j]);
                  N_tensor_array[6*i+5] +=  factor*(dipole::internal::rij_tensor_zz[lc][j]);

               }
            }
         }
         // // Print tensor: uncomment if you want to check the tensor components
         // std::cout << "*----------------------------------*" << std::endl;
         // std::cout << "lc = " << lc << "\ti = " << i << "\tNat_cell_i = " << dipole::internal::cells_num_atoms_in_cell[i]  << "\tN_self_i = " << N_self_array[i] << std::endl;
         // std::cout << N_tensor_array[6*i+0] << "\t" << N_tensor_array[6*i+1] << "\t" << N_tensor_array[6*i+2] << "\n";
         // std::cout << N_tensor_array[6*i+1] << "\t" << N_tensor_array[6*i+3] << "\t" << N_tensor_array[6*i+4] << "\n";
         // std::cout << N_tensor_array[6*i+2] << "\t" << N_tensor_array[6*i+4] << "\t" << N_tensor_array[6*i+5] << "\n";
         // std::cout << "*----------------------------------*" << std::endl;
         // std::cout << std::endl;
      }

      // if vampire is running in parallel, all cpus send demag field to root proc
      #ifdef MPICF
         dipole::internal::send_cells_demag_factor(cells::cell_id_array,
                                                   N_tensor_array,
                                                   cells::num_local_cells
         );
         // Could potentially replace dipole::internal::send_cells_demag_factor() with MPI_Reduce and MPI_MAX op, but depends on implementation.
         // Leaving as-is for now
      #endif

      // Compute demag factor only on root process
      if (vmpi::my_rank == 0)
      {

         int cells_non_zero = 0;          /// Counter for cells with atoms inside
         int num_atoms_magnetic = 0.0;   // Initialise tot num of magnetic atoms
         double Vtot = 0.0;               /// Initialise total volume of the system
         // Define and initialise Demag factor N tensor components
         double Nxx = 0.0;
         double Nxy = 0.0;
         double Nxz = 0.0;
         double Nyy = 0.0;
         double Nyz = 0.0;
         double Nzz = 0.0;

         // Calculate number of magnetic atoms and Vtot
         for(int i = 0; i < cells::num_cells; i++)
         {
            if (dipole::internal::cells_num_atoms_in_cell[i] > 0)
            {
               num_atoms_magnetic += dipole::internal::cells_num_atoms_in_cell[i];
               cells_non_zero ++;
            }
         }

      //   for(int lc=0; lc<dipole::internal::cells_num_local_cells; lc++){
      //      int i = cells::cell_id_array[lc];
         for(int i=0; i<cells::num_cells; i++){
            if(dipole::internal::cells_num_atoms_in_cell[i]>0){
               // Calculate total volume summing over cells volume
               Vtot += dipole::internal::cells_volume_array[i]; // /double(dipole::internal::cells_num_atoms_in_cell[i]); //+= dipole::internal::cells_volume_array[i];
               // Calculate Demag tensor
               Nxx += N_tensor_array[6*i+0];
               Nxy += N_tensor_array[6*i+1];
               Nxz += N_tensor_array[6*i+2];
               Nyy += N_tensor_array[6*i+3];
               Nyz += N_tensor_array[6*i+4];
               Nzz += N_tensor_array[6*i+5];
               // // Print tensor: uncomment if you want to check the tensor components
               // std::cout << "*----------- total tensor -------------*" << std::endl;
               // std::cout << "i = " << i << "\tNat_cell_i = " << dipole::internal::cells_num_atoms_in_cell[i]  << "\tN_self_i = " << N_self_array[i] << std::endl;
               // //std::cout << "lc = " << lc << "\ti = " << i << "\tNat_cell_i = " << dipole::internal::cells_num_atoms_in_cell[i]  << "\tN_self_i = " << N_self_array[i] << std::endl;
               // std::cout << Nxx << "\t" << Nxy << "\t" << Nxz << "\n";
               // std::cout << Nxy << "\t" << Nyy << "\t" << Nyz << "\n";
               // std::cout << Nxz << "\t" << Nyz << "\t" << Nzz << "\n";
               // std::cout << "*----------------------------------*" << std::endl;
               // std::cout << std::endl;
            }
         }

         // Normalise N by the total number of magnetic atoms ~ volume
         Nxx =  Nxx / ( double(num_atoms_magnetic) * double(num_atoms_magnetic) );
         Nxy =  Nxy / ( double(num_atoms_magnetic) * double(num_atoms_magnetic) );
         Nxz =  Nxz / ( double(num_atoms_magnetic) * double(num_atoms_magnetic) );
         Nyy =  Nyy / ( double(num_atoms_magnetic) * double(num_atoms_magnetic) );
         Nyz =  Nyz / ( double(num_atoms_magnetic) * double(num_atoms_magnetic) );
         Nzz =  Nzz / ( double(num_atoms_magnetic) * double(num_atoms_magnetic) );

         // // Print tensor: uncomment if you want to check the tensor components
         // std::cout << "\n*----------- total tensor -------------*" << std::endl;
         // // std::cout << "i = " << i << "\tNat_cell_i = " << dipole::internal::cells_num_atoms_in_cell[i]  << "\tN_self_i = " << N_self_array[i] << std::endl;
         // // //std::cout << "lc = " << lc << "\ti = " << i << "\tNat_cell_i = " << dipole::internal::cells_num_atoms_in_cell[i]  << "\tN_self_i = " << N_self_array[i] << std::endl;
         // std::cout << "\ncells_non_zero\t" << cells_non_zero << std::endl;
         // std::cout << "\nVtot\t" << Vtot << "\tVtot/Nmag\t" << Vtot/double(num_atoms_magnetic) << std::endl;
         // std::cout << Nxx << "\t" << Nxy << "\t" << Nxz << "\n";
         // std::cout << Nxy << "\t" << Nyy << "\t" << Nyz << "\n";
         // std::cout << Nxz << "\t" << Nyz << "\t" << Nzz << "\n";
         // std::cout << "*----------------------------------*" << std::endl;
         // std::cout << std::endl;

         Nxx =  -(Nxx * Vtot)/(4.0*M_PI) + 0.3333333333333;
         Nxy =  -(Nxy * Vtot)/(4.0*M_PI) ;
         Nxz =  -(Nxz * Vtot)/(4.0*M_PI) ;
         Nyy =  -(Nyy * Vtot)/(4.0*M_PI) + 0.3333333333333;
         Nyz =  -(Nyz * Vtot)/(4.0*M_PI) ;
         Nzz =  -(Nzz * Vtot)/(4.0*M_PI) + 0.3333333333333;

         // Write demag factor to log file zlog << zTs() <<
         zlog << zTs() << "Demagnetisation tensor in format Nxx   Nxy   Nxz   Nyx   Nyy   Nyz   Nzx   Nzy   Nzz :\n";
         zlog << zTs() << Nxx << "\t" << Nxy << "\t" << Nxz << "\t";
         zlog <<          Nxy << "\t" << Nyy << "\t" << Nyz << "\t";
         zlog <<          Nxz << "\t" << Nyz << "\t" << Nzz << "\n";
      }
      // Clear memory
      N_tensor_array.clear();

	   //--------------------------------------------------/
      // End of output dipolar tensor
	   //--------------------------------------------------/

    	return;
    }
} // end of dipole namespace
