//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Rory Pond and Richard F L Evans 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
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
#include "constants.hpp"
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

/// @brief Atomistic output function
///
/// @details Outputs formatted data snapshot for visualisation
///
///   #------------------------------------------------------
///   # Atomistic spin configuration file for vampire
///   #------------------------------------------------------
///   # Date: xx/xx/xxxx xx.xx.xx
///   #------------------------------------------------------
///   Number of spins: $n_spins
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
///   Number of spin files: $n_files
///   atoms-000ID000-00CPU0.cfg
///   atoms-000ID000-00CPU1.cfg
///   atoms-000ID000-00CPU2.cfg
///   #------------------------------------------------------
///   Number of local spins: $n_loc_spins
///   $sx $sy $sz
///   $sx $sy ...
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
///   Created: 30/05/2011
///   Revision:   ---
///=====================================================================================
///
double legacy_atoms(){

   // instantiate timer
   vutil::vtimer_t timer;

   // check calling of routine if error checking is activated
   if (err::check == true)
   {
      std::cout << "config::atoms has been called" << std::endl;
   }

   /* #ifdef MPICF
      const int num_atoms = vmpi::num_core_atoms + vmpi::num_bdry_atoms;
   #else
      const int num_atoms = atoms::num_atoms;
   #endif */ //unused variable

   // Set local output filename
   std::stringstream file_sstr;
   file_sstr << "atoms-";
   // Set CPUID on non-root process
   if (vmpi::my_rank != 0)
   {
      file_sstr << std::setfill('0') << std::setw(5) << vmpi::my_rank << "-";
   }
   file_sstr << std::setfill('0') << std::setw(8) << sim::output_atoms_file_counter;
   file_sstr << ".cfg";
   std::string cfg_file = file_sstr.str();
   const char *cfg_filec = cfg_file.c_str();

   // Declare and open output file
   std::ofstream cfg_file_ofstr;
   cfg_file_ofstr.open(cfg_filec);

   // Output masterfile header on root process
   if (vmpi::my_rank == 0)
   {
      // Get system date
      time_t rawtime = time(NULL);
      struct tm *timeinfo = localtime(&rawtime);

      cfg_file_ofstr << "#------------------------------------------------------" << std::endl;
      cfg_file_ofstr << "# Atomistic spin configuration file for vampire" << std::endl;
      cfg_file_ofstr << "#------------------------------------------------------" << std::endl;
      cfg_file_ofstr << "# Date: " << asctime(timeinfo);
      cfg_file_ofstr << "#------------------------------------------------------" << std::endl;
      cfg_file_ofstr << "Number of spins: " << config::internal::total_output_atoms << std::endl;
      cfg_file_ofstr << "System dimensions:" << cs::system_dimensions[0] << "\t" << cs::system_dimensions[1] << "\t" << cs::system_dimensions[2] << std::endl;
      cfg_file_ofstr << "Coordinates-file: atoms-coord.cfg" << std::endl;
      cfg_file_ofstr << "Time: " << double(sim::time) * mp::dt_SI << std::endl;
      cfg_file_ofstr << "Field: " << sim::H_applied << std::endl;
      cfg_file_ofstr << "Temperature: " << sim::temperature << std::endl;
      cfg_file_ofstr << "Magnetisation: " << stats::system_magnetization.output_normalized_magnetization(false) << std::endl;
      cfg_file_ofstr << "Number of Materials: " << mp::num_materials << std::endl;
      for (int mat = 0; mat < mp::num_materials; mat++)
      {
         cfg_file_ofstr << mp::material[mat].mu_s_SI << std::endl;
      }
      cfg_file_ofstr << "#------------------------------------------------------" << std::endl;
      cfg_file_ofstr << "Number of spin files: " << vmpi::num_processors - 1 << std::endl;
      for (int p = 1; p < vmpi::num_processors; p++)
      {
         std::stringstream cfg_sstr;
         cfg_sstr << "atoms-" << std::setfill('0') << std::setw(5) << p << "-" << std::setfill('0') << std::setw(8) << sim::output_atoms_file_counter << ".cfg";
         cfg_file_ofstr << cfg_sstr.str() << std::endl;
      }
      cfg_file_ofstr << "#------------------------------------------------------" << std::endl;
   }

   // start timer
   timer.start();

   // Everyone now outputs their atom list
   cfg_file_ofstr << local_output_atom_list.size() << std::endl;
   for (size_t i = 0; i < local_output_atom_list.size(); i++)
   {
      const int atom = local_output_atom_list[i];
      cfg_file_ofstr << atoms::x_spin_array[atom] << "\t" << atoms::y_spin_array[atom] << "\t" << atoms::z_spin_array[atom] << std::endl;
   }

   // stop timer
   timer.stop();

   // close output file
   cfg_file_ofstr.close();

   double io_time = timer.elapsed_time(); // seconds

   #ifdef MPICF
      // find maximum time for i/o
      double max_io_time = 0.0;
      // calculate actual bandwidth on root process
      MPI_Reduce(&io_time, &max_io_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      io_time = max_io_time;
   #endif

   return io_time;

}

/// @brief Atomistic output function
///
/// @details Outputs formatted data snapshot for visualisation
///
///   //------------------------------------------------------
///   // Atomistic coordinate configuration file for vampire
///   //------------------------------------------------------
///   // Date: xx/xx/xxxx xx.xx.xx
///   //------------------------------------------------------
///   Number of atoms: $n_spins
///   //------------------------------------------------------
///   Number of atom files: $n_files
///   atoms-coords-00CPU0.cfg
///   atoms-coords-00CPU1.cfg
///   atoms-coords-00CPU2.cfg
///   //------------------------------------------------------
///   Number of local atoms: $n_loc_atoms
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
///   Created: 31/05/2011
///   Revision:   ---
///=====================================================================================
///
double legacy_atoms_coords()
{

   // instantiate timer
   vutil::vtimer_t timer;

   // check calling of routine if error checking is activated
   if (err::check == true)
   {
      std::cout << "config::atoms_coords has been called" << std::endl;
   }

   // Set local output filename
   std::stringstream file_sstr;
   file_sstr << "atoms-coords";
   // Set CPUID on non-root process
   if (vmpi::my_rank != 0)
   {
      file_sstr << "-" << std::setfill('0') << std::setw(5) << vmpi::my_rank;
   }
   file_sstr << ".cfg";
   std::string cfg_file = file_sstr.str();
   const char *cfg_filec = cfg_file.c_str();

   // Declare and open output file
   std::ofstream cfg_file_ofstr;
   cfg_file_ofstr.open(cfg_filec);

   // start timer
   timer.start();

   // Output masterfile header on root process
   if (vmpi::my_rank == 0)
   {
      // Get system date
      time_t rawtime = time(NULL);
      struct tm *timeinfo = localtime(&rawtime);

      cfg_file_ofstr << "#------------------------------------------------------" << std::endl;
      cfg_file_ofstr << "# Atomistic coordinates configuration file for vampire" << std::endl;
      cfg_file_ofstr << "#------------------------------------------------------" << std::endl;
      cfg_file_ofstr << "# Date: " << asctime(timeinfo);
      cfg_file_ofstr << "#------------------------------------------------------" << std::endl;
      cfg_file_ofstr << "Number of atoms: " << config::internal::total_output_atoms << std::endl;
      cfg_file_ofstr << "#------------------------------------------------------" << std::endl;
      cfg_file_ofstr << "Number of spin files: " << vmpi::num_processors - 1 << std::endl;
      for (int p = 1; p < vmpi::num_processors; p++)
      {
         std::stringstream cfg_sstr;
         cfg_sstr << "atoms-coords-" << std::setfill('0') << std::setw(5) << p << ".cfg";
         cfg_file_ofstr << cfg_sstr.str() << std::endl;
      }
      cfg_file_ofstr << "#------------------------------------------------------" << std::endl;
   }

   // Everyone now outputs their atom list
   cfg_file_ofstr << local_output_atom_list.size() << std::endl;
   for (size_t i = 0; i < local_output_atom_list.size(); i++)
   {
      const int atom = local_output_atom_list[i];
      cfg_file_ofstr << atoms::type_array[atom] << "\t" << atoms::category_array[atom] << "\t" << atoms::x_coord_array[atom] << "\t" << atoms::y_coord_array[atom] << "\t" << atoms::z_coord_array[atom] << "\t";
      if (config::internal::identify_surface_atoms == true && atoms::surface_array[atom] == true)
         cfg_file_ofstr << "O " << std::endl;
      else
         cfg_file_ofstr << mp::material[atoms::type_array[atom]].element << std::endl;
   }

   // stop the timer
   timer.stop();

   cfg_file_ofstr.close();

   double io_time = timer.elapsed_time(); // seconds

   #ifdef MPICF
      // find maximum time for i/o
      double max_io_time = 0.0;
      // calculate actual bandwidth on root process
      MPI_Reduce(&io_time, &max_io_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      io_time = max_io_time;
   #endif

   return io_time;

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
void legacy_cells()
{

   // wait for all processes
   vmpi::barrier();

   // update cells magnetization
   cells::mag();

   // instantiate timer
   vutil::vtimer_t timer;

   // Set local output filename
   std::stringstream file_sstr;
   file_sstr << "cells-";
   file_sstr << std::setfill('0') << std::setw(8) << sim::output_cells_file_counter;
   file_sstr << ".cfg";
   std::string cfg_file = file_sstr.str();
   const char *cfg_filec = cfg_file.c_str();

   #ifdef MPICF
      // if flag to print cells field is active, all cpus send cells field to root proc
      if(dipole::activated) dipole::send_cells_field(cells::cell_id_array,
                                                     dipole::cells_field_array_x,
                                                     dipole::cells_field_array_y,
                                                     dipole::cells_field_array_z,
                                                     cells::num_local_cells);
   #endif

   // start timer
   timer.start();

   // Output masterfile header on root process
   if (vmpi::my_rank == 0)
   {

      zlog << zTs() << "Outputting cell configuration " << sim::output_cells_file_counter << " to disk" << std::flush;

      // Declare and open output file
      std::ofstream cfg_file_ofstr;
      cfg_file_ofstr.open(cfg_filec);

      // Get system date
      time_t rawtime = time(NULL);
      struct tm *timeinfo = localtime(&rawtime);

      cfg_file_ofstr << "#------------------------------------------------------" << std::endl;
      cfg_file_ofstr << "# Cell configuration file for vampire" << std::endl;
      cfg_file_ofstr << "#------------------------------------------------------" << std::endl;
      cfg_file_ofstr << "# Date: " << asctime(timeinfo);
      cfg_file_ofstr << "#------------------------------------------------------" << std::endl;
      cfg_file_ofstr << "# Number of spins: " << cells::num_cells << std::endl;
      cfg_file_ofstr << "# System dimensions:" << cs::system_dimensions[0] << "\t" << cs::system_dimensions[1] << "\t" << cs::system_dimensions[2] << std::endl;
      cfg_file_ofstr << "# Coordinates-file: cells-coord.cfg" << std::endl;
      cfg_file_ofstr << "# Time: " << double(sim::time) * mp::dt_SI << std::endl;
      cfg_file_ofstr << "# Field: " << sim::H_applied << std::endl;
      cfg_file_ofstr << "# Temperature: " << sim::temperature << std::endl;
      cfg_file_ofstr << "# Magnetisation: " << stats::system_magnetization.output_normalized_magnetization(false) << std::endl;
      cfg_file_ofstr << "#------------------------------------------------------" << std::endl;

      const double inv_muB = 1.0 / constants::muB;

      // Root process now outputs the cell magnetisations
      for (int cell = 0; cell < cells::num_cells; cell++){
         // get cells magnetization
         const double mx = cells::mag_array_x[cell]*inv_muB;
         const double my = cells::mag_array_y[cell]*inv_muB;
         const double mz = cells::mag_array_z[cell]*inv_muB;
         const double mm = sqrt(mx*mx + my*my + mz*mz); // actual vector length (Bohr magnetons)
         const double imm = 1.0/mm;
         const double rm = mm / cells::pos_and_mom_array[4*cell+3]/inv_muB; // relative moment

         // only output cells with magnetic moments
         if (cells::num_atoms_in_cell_global[cell] > 0){
            cfg_file_ofstr << mx*imm << "\t" << my*imm << "\t" << mz*imm << "\t" << rm << "\t" << mm << "\t";
            if(dipole::activated) cfg_file_ofstr << dipole::cells_field_array_x[cell] << "\t" << dipole::cells_field_array_y[cell] << "\t" << dipole::cells_field_array_z[cell] << "\n";
            else cfg_file_ofstr << "\n";
         }

      }

      cfg_file_ofstr.close();
   }

   // stop the timer
   timer.stop();

   double data_size = double(config::internal::total_output_cells) * 3.0 * sizeof(double);
   if(dipole::activated) data_size = data_size * 2.0;
   const double io_time = timer.elapsed_time();

   zlog << " of size " << data_size*1.0e-6 << " MB [ " << data_size*1.0e-9/timer.elapsed_time() << " GB/s in " << io_time << " s ]" << std::endl;

   sim::output_cells_file_counter++;

   // wait for all processes
   vmpi::barrier();

   return;
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
void legacy_cells_coords()
{

   // wait for all processes
   vmpi::barrier();

   // instantiate timer
   vutil::vtimer_t timer;

   // Set local output filename
   std::stringstream file_sstr;
   file_sstr << "cells-coords";
   file_sstr << ".cfg";
   std::string cfg_file = file_sstr.str();
   const char *cfg_filec = cfg_file.c_str();

   // start timer
   timer.start();

   // Output masterfile header on root process
   if (vmpi::my_rank == 0)
   {

      zlog << zTs() << "Outputting cell coordinates to disk" << std::flush;

      // Declare and open output file
      std::ofstream cfg_file_ofstr;
      cfg_file_ofstr.open(cfg_filec);

      // Get system date
      time_t rawtime = time(NULL);
      struct tm *timeinfo = localtime(&rawtime);

      cfg_file_ofstr << "#------------------------------------------------------" << std::endl;
      cfg_file_ofstr << "# Cell coordinates configuration file for vampire" << std::endl;
      cfg_file_ofstr << "#------------------------------------------------------" << std::endl;
      cfg_file_ofstr << "# Date: " << asctime(timeinfo);
      cfg_file_ofstr << "#------------------------------------------------------" << std::endl;
      cfg_file_ofstr << "# Number of cells: " << cells::num_cells << std::endl;
      cfg_file_ofstr << "#------------------------------------------------------" << std::endl;
      cfg_file_ofstr << "# cell size: " << cells::macro_cell_size_x << "\t" << cells::macro_cell_size_y << "\t" << cells::macro_cell_size_z << "\t" <<std::endl;
      cfg_file_ofstr << "#------------------------------------------------------" << std::endl;
      cfg_file_ofstr << "#" << std::endl;
      cfg_file_ofstr << "#" << std::endl;
      cfg_file_ofstr << "#" << std::endl;
      cfg_file_ofstr << "#------------------------------------------------------" << std::endl;

      config::internal::total_output_cells = 0;

      const double inv_muB = 1.0 / constants::muB;

      for (int cell = 0; cell < cells::num_cells; cell++){
         // only output cells with magnetic moments
         if (cells::num_atoms_in_cell_global[cell] > 0){
            const double nm = cells::num_atoms_in_cell_global[cell];
            const double cx = cells::pos_and_mom_array[4 * cell + 0];
            const double cy = cells::pos_and_mom_array[4 * cell + 1];
            const double cz = cells::pos_and_mom_array[4 * cell + 2];
            const double mm = cells::pos_and_mom_array[4 * cell + 3]*inv_muB;
            // calculate cell corners
            const double cxp = cx + cells::macro_cell_size_x * 0.5;
            const double cyp = cy + cells::macro_cell_size_y * 0.5;
            const double czp = cz + cells::macro_cell_size_z * 0.5;
            const double cxm = cx - cells::macro_cell_size_x * 0.5;
            const double cym = cy - cells::macro_cell_size_y * 0.5;
            const double czm = cz - cells::macro_cell_size_z * 0.5;

            cfg_file_ofstr << cell << "\t" << nm << "\t" << mm  << "\t"
                           << cx  << "\t" << cy  << "\t" << cz  << "\t"
                           << cxm << "\t" << cym << "\t" << czm << "\t"
                           << cxp << "\t" << cyp << "\t" << czp << std::endl;

            // increment cell counter
            config::internal::total_output_cells++;

         }
      }
      cfg_file_ofstr.close();
   }

   // stop the timer
   timer.stop();

   const double data_size = double(config::internal::total_output_cells) * 2.0 * sizeof(int) * 3.0 * sizeof(double);
   const double io_time = timer.elapsed_time();

   zlog << " of size " << data_size*1.0e-6 << " MB [ " << data_size*1.0e-9/timer.elapsed_time() << " GB/s in " << io_time << " s ]" << std::endl;

   // wait for all processes
   vmpi::barrier();

   return;
}
}
}
