//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) rory.pond 2016. All rights reserved.
//
//   Email: rory.pond@york.ac.uk
//
//------------------------------------------------------------------------------
//

#include <stdio.h>

// C++ standard library headers
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

// Vampire headers
#include "atoms.hpp"
#include "cells.hpp"
#include "dipole.hpp"
#include "errors.hpp"
#include "LLG.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "vutil.hpp"

// config module headers
#include "config.hpp"
#include "internal.hpp"

namespace config
{

namespace internal
{
std::vector<float> total_buffer;
std::vector<float> localbuffer;
/// @brief Atomistic output function
///
/// @details Outputs formatted data snapshot for visualisation
///
///	#------------------------------------------------------
///	# Atomistic spin configuration file for vampire
///	#------------------------------------------------------
///	# Date: xx/xx/xxxx xx.xx.xx
///	#------------------------------------------------------
///	Number of spins: $n_spins
///	System dimensions: $max_x $max_y $max_z
///	Coordinates-file: $coord_file
///	Time: $t
///	Field: $H
///	Temperature: $T
///	Magnetisation: $mx $my $mz
///	Number of Materials: $n_mat
///	Material Properties 1:	$mu_s	$mmx $mmy $mmz $mm
///	Material Properties 2:	$mu_s	$mmx $mmy $mmz ...
///	#------------------------------------------------------
///	Number of spin files: $n_files
///	atoms-000ID000-00CPU0.cfg
///	atoms-000ID000-00CPU1.cfg
///	atoms-000ID000-00CPU2.cfg
///	#------------------------------------------------------
///	Number of local spins: $n_loc_spins
///	$sx $sy $sz
///	$sx $sy ...
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    30/05/2011
///
/// @internal
///	Created:		30/05/2011
///	Revision:	---
///=====================================================================================
///

void atoms_new()
{

   // instantiate timer
   vutil::vtimer_t timer;

   // start timer
   timer.start();

   // check calling of routine if error checking is activated
   if (err::check == true)
   {
      std::cout << "config::atoms has been called" << std::endl;
   }

   // Output Meta Data
      const std::vector<double> m_l = stats::system_magnetization.get_magnetization();

      int files;

      if(output_all)
         files = vmpi::num_processors;

      if(output_gather)
         files = vmpi::num_io_processors;

      if (output_mpi_io)
         files = 1;
      

      if (vmpi::my_rank == 0)
         write_meta( double(sim::time) * mp::dt_SI , // time (seconds)
                     sim::temperature, // system temperature (Kelvin)
                     sim::H_vec[0], // applied field components (Tesla)
                     sim::H_vec[1],
                     sim::H_vec[2],
                     m_l[0], // magnetization components (normalized)
                     m_l[0],
                     m_l[0],
                     files); // number of files

      //Output spindata
      int local_size = local_output_atom_list.size() * 3;
      localbuffer.clear();
      total_buffer.clear();
      localbuffer.resize(local_size);
      copy_data_to_buffer(atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array, local_output_atom_list, localbuffer);

   #ifdef MPICF
      if (output_mpi_io)
      {
            MPI_File fh;
            MPI_Status status; 
            std::string filestring = data_filename(false);
            char *filename = (char*)filestring.c_str();
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
            MPI_File_write_ordered(fh, &localbuffer[0], local_size, MPI_FLOAT, &status);
            MPI_File_close(&fh);
      }else if(output_all)
      {
         write_data(localbuffer, false);
      }else if(output_gather)
      {
            
            int total_size=0;
            MPI_Reduce(&local_size, &total_size, 1, MPI_INT, MPI_SUM, vmpi::io_processor,vmpi::io_comm);
            total_buffer.resize(total_size);
            //MPI_Gather( &localbuffer[0], local_size, MPI_FLOAT, &total_buffer[0], total_size, MPI_FLOAT, vmpi::io_processor, vmpi::io_comm);
            
            std::vector<int> buffer_sizes(vmpi::size_io_group), buffer_offset(vmpi::size_io_group);
            MPI_Gather( &local_size, 1, MPI_INT, &buffer_sizes[0], vmpi::size_io_group, MPI_INT, vmpi::io_processor, vmpi::io_comm);
            int temp=0;
            for (int i = 0; i < vmpi::size_io_group; i++)
            {
               buffer_offset[i]=temp;
               temp+=buffer_sizes[i];
               total_size += buffer_sizes[i];
            }
             
            //MPI_Gather( &localbuffer[0], local_size, MPI_FLOAT, &total_buffer[0], total_size, MPI_FLOAT, vmpi::io_processor, vmpi::io_comm);
            //MPI_Gatherv( &localbuffer[0], local_size, MPI_FLOAT, &total_buffer[0], &buffer_sizes[0], &buffer_offset[0], MPI_FLOAT, vmpi::io_processor, vmpi::io_comm);
            
            if (vmpi::my_io_rank == vmpi::io_processor)
               write_data(total_buffer, false);
      }
   #else
      write_data(localbuffer, false);
   #endif

      // stop the timer
   double local_time = timer.elapsed_time(); // seconds

      // get file size (bytes)
   double local_data_size = double(sizeof(float) * local_size);
#ifdef MPICF
// aggregate bandwidth
   
   double total_time;
   MPI_Reduce(&local_time, &total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
   double total_data_size;
   MPI_Reduce(&local_data_size, &total_data_size, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if(vmpi::my_rank==0)
      zlog << "Data transfer time " << total_time << std::endl;
      zlog << "Data transfer speed " << 1.0e-6 * total_data_size / total_time << " MB/s" << std::endl;
#else
   // calculate data rate and output to log
   zlog << "Data transfer time "  << local_time << std::endl;
   zlog << "Data transfer speed " << 1.0e-6 * local_data_size / local_time << " MB/s" << std::endl;
#endif
}

/// @brief Atomistic output function
///
/// @details Outputs formatted data snapshot for visualisation
///
///	//------------------------------------------------------
///	// Atomistic coordinate configuration file for vampire
///	//------------------------------------------------------
///	// Date: xx/xx/xxxx xx.xx.xx
///	//------------------------------------------------------
///	Number of atoms: $n_spins
///	//------------------------------------------------------
///	Number of atom files: $n_files
///	atoms-coords-00CPU0.cfg
///	atoms-coords-00CPU1.cfg
///	atoms-coords-00CPU2.cfg
///	//------------------------------------------------------
///	Number of local atoms: $n_loc_atoms
///	$material $category $x $y $z $species
///	$material $category $x $y $z $species
///	$material $category $x ...
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    31/05/2011
///
/// @internal
///	Created:		31/05/2011
///	Revision:	---
///=====================================================================================
///
void atoms_coords_new()
{
      // instantiate timer
   vutil::vtimer_t timer;

   // start timer
   timer.start();

   // check calling of routine if error checking is activated
   if (err::check == true)
   {
      std::cout << "config::atoms_coords has been called" << std::endl;
   }

   #ifdef MPICF
      const int num_atoms = vmpi::num_core_atoms + vmpi::num_bdry_atoms;
   #else
      const int num_atoms = atoms::num_atoms;
   #endif

   // resize atom list to zero
   local_output_atom_list.resize(0);

   // get output bounds
   double minB[3] = {atoms_output_min[0] * cs::system_dimensions[0],
                     atoms_output_min[1] * cs::system_dimensions[1],
                     atoms_output_min[2] * cs::system_dimensions[2]};

   double maxB[3] = {atoms_output_max[0] * cs::system_dimensions[0],
                     atoms_output_max[1] * cs::system_dimensions[1],
                     atoms_output_max[2] * cs::system_dimensions[2]};

   // loop over all local atoms and record output list
   for (int atom = 0; atom < num_atoms; atom++)
   {

      const double cc[3] = {atoms::x_coord_array[atom], atoms::y_coord_array[atom], atoms::z_coord_array[atom]};

      // check atom within output bounds
      if ((cc[0] >= minB[0]) && (cc[0] <= maxB[0]))
      {
         if ((cc[1] >= minB[1]) && (cc[1] <= maxB[1]))
         {
            if ((cc[2] >= minB[2]) && (cc[2] <= maxB[2]))
            {
               local_output_atom_list.push_back(atom);
            }
         }
      }
   }

   // calculate total atoms to output
   #ifdef MPICF
      int local_atoms = local_output_atom_list.size();
      int total_atoms;
      //std::cerr << vmpi::my_rank << "\t" << local_atoms << &local_atoms << "\t" << &total_atoms << std::endl;
      //MPI::COMM_WORLD.Barrier();
      MPI::COMM_WORLD.Allreduce(&local_atoms, &total_atoms, 1, MPI_INT, MPI_SUM);
      total_output_atoms = total_atoms;
   //std::cerr << vmpi::my_rank << "\t" << total_atoms << "\t" << &local_atoms << "\t" << &total_atoms << std::endl;
   //MPI::COMM_WORLD.Barrier();
   #else
      total_output_atoms = local_output_atom_list.size();
   #endif

   // Output Meta Data
   write_coordinate_meta();

      //Output Coord Data
      int local_size = local_output_atom_list.size() * 3;
      localbuffer.clear();
      total_buffer.clear();
      localbuffer.resize(local_size);
      copy_data_to_buffer(atoms::x_coord_array, atoms::y_coord_array, atoms::z_coord_array, local_output_atom_list, localbuffer);

   #ifdef MPICF
      if (output_mpi_io)
      {
            MPI_File fh;
            MPI_Status status;
            std::string filestring = data_filename(true);
            char *filename = (char*)filestring.c_str();
            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
            MPI_File_write_ordered(fh, &localbuffer[0], local_size, MPI_FLOAT, &status);
            MPI_File_close(&fh);
      }else if(output_all)
      {
         write_data(localbuffer, true);
      }else if(output_gather)
      {
            /*printf("test 1 rank %d", vmpi::my_rank);
            int total_size=0;
            MPI_Reduce(&local_size, &total_size, 1, MPI_INT, MPI_SUM, vmpi::io_processor,vmpi::io_comm);
            MPI_Barrier( vmpi::io_comm );
            total_buffer.resize(total_size);
            printf("The IO Processor is %d \n", vmpi::io_processor);
            //std::vector<int> buffer_sizes(vmpi::size_io_group), buffer_offset(vmpi::size_io_group);
            int buffer_sizes[vmpi::size_io_group];
            int buffer_offset[vmpi::size_io_group];
            MPI_Gather( &local_size, 1, MPI_INT, &buffer_sizes[0], vmpi::size_io_group, MPI_INT, vmpi::io_processor, MPI_COMM_WORLD);//vmpi::io_comm);
            int temp=0;
            for (int i = 0; i < vmpi::size_io_group; i++)
            {
               buffer_offset[i]=temp;
               temp+=buffer_sizes[i];
            }
               printf("\n rank %d, local size %d", vmpi::my_rank, local_size);
            if (vmpi::my_io_rank == vmpi::io_processor)
               printf("\n local [0] = %d, local [1] = %d, local [2] = %d, local [3] = %d \n" , buffer_sizes[0], buffer_sizes[1], buffer_sizes[2], buffer_sizes[3]);
            //MPI_Gather( &localbuffer[0], local_size, MPI_FLOAT, &total_buffer[0], total_size, MPI_FLOAT, vmpi::io_processor, vmpi::io_comm);
            //MPI_Gatherv( &localbuffer[0], local_size, MPI_FLOAT, &total_buffer[0], &buffer_sizes[0], &buffer_offset[0], MPI_FLOAT, vmpi::io_processor, vmpi::io_comm);
            
            if (vmpi::my_io_rank == vmpi::io_processor)
               write_data(total_buffer, true);
            */
      }
   #else
      write_data(localbuffer, true);

   #endif

         // stop the timer
   double local_time = timer.elapsed_time(); // seconds

      // get file size (bytes)
   double local_data_size = double(sizeof(float) * local_size);
#ifdef MPICF
// aggregate bandwidth
   
   double total_time;
   MPI_Reduce(&local_time, &total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
   double total_data_size;
   MPI_Reduce(&local_data_size, &total_data_size, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if(vmpi::my_rank==0)
      zlog << "Data transfer time " << total_time << std::endl;
      zlog << "Data transfer speed " << 1.0e-6 * total_data_size / total_time << " MB/s" << std::endl;
#else
   // calculate data rate and output to log
   zlog << "Data transfer time "  << local_time << std::endl;
   zlog << "Data transfer speed " << 1.0e-6 * local_data_size / local_time << " MB/s" << std::endl;
#endif
}


/// @brief Cell output function
///
/// @details Outputs formatted data snapshot for visualisation
///
///   #------------------------------------------------------
///   # Cell configuration file for vampire
///   #------------------------------------------------------
///   # Date: xx/xx/xxxx xx.xx.xx
///   #------------------------------------------------------
///   Number of cells: $n_cells
///   System dimensions: $max_x $max_y $max_z
///   Coordinates-file: $coord_file
///   Time: $t
///   Field: $H
///   Temperature: $T
///   Magnetisation: $mx $my $mz
///   Number of Materials: $n_mat
///   Material Properties 1:  $mu_s $mmx $mmy $mmz $mm
///   Material Properties 2:  $mu_s $mmx $mmy $mmz ...
///   #------------------------------------------------------
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    26/04/2013
///
/// @internal
///   Created:    26/04/2013
///   Revision:     ---
///=====================================================================================
///
void cells_new()
{

}

/// @brief Cells output function
///
/// @details Outputs formatted data snapshot for visualisation
///
///   //------------------------------------------------------
///   // Atomistic coordinate configuration file for vampire
///   //------------------------------------------------------
///   // Date: xx/xx/xxxx xx.xx.xx
///   //------------------------------------------------------
///   Number of cells: $n_cells
///   //------------------------------------------------------
///   Number of atom files: $n_files
///   atoms-coords-00CPU0.cfg
///   atoms-coords-00CPU1.cfg
///   atoms-coords-00CPU2.cfg
///   //------------------------------------------------------
///   Number of local cells: $n_loc_cells
///   $material $category $x $y $z $species
///   $material $category $x $y $z $species
///   $material $category $x ...
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    31/05/2011
///
/// @internal
///   Created:    31/05/2011
///   Revision:     ---
///=====================================================================================
///
void cells_coords_new()
{
      
}
}
}