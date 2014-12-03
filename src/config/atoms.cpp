//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <sstream>

// Vampire headers
#include "errors.hpp"
#include "ltmp.hpp"
#include "vio.hpp"

// Localised temperature pulse headers
#include "internal.hpp"

// calculate pitch - number of processors per output node

// last processor in batch given job of managing data output

// set up MPI routines and storage

// cast spin array to float

// all group processors send data to my_io_rank in non-blocking fashion

// next time routine is called write data to disk (remember to code for check written function)

// data blocked in [sx][sy][sz]

// file format
// n_spins
// [sx][sy][sz]











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
///   Created:    30/05/2011
///   Revision:     ---
///=====================================================================================
///
   void atoms(){

      // check calling of routine if error checking is activated
      if(err::check==true){std::cout << "vout::atoms has been called" << std::endl;}

      #ifdef MPICF
         const int num_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
      #else
         const int num_atoms = atoms::num_atoms;
      #endif

      // Set local output filename
      std::stringstream file_sstr;
      file_sstr << "atoms-";
      // Set CPUID on non-root process
      if(vmpi::my_rank!=0){
         file_sstr << std::setfill('0') << std::setw(5) << vmpi::my_rank << "-";
      }
      file_sstr << std::setfill('0') << std::setw(8) << output_atoms_file_counter;
      file_sstr << ".cfg";
      std::string cfg_file = file_sstr.str();
      const char* cfg_filec = cfg_file.c_str();

      // Output informative message to log file
      zlog << zTs() << "Outputting configuration file " << cfg_file << " to disk" << std::endl;

      // Declare and open output file
      std::ofstream cfg_file_ofstr;
      cfg_file_ofstr.open (cfg_filec);

      // Output masterfile header on root process
      if(vmpi::my_rank==0){
         // Get system date
      time_t rawtime = time(NULL);
      struct tm * timeinfo = localtime(&rawtime);

         cfg_file_ofstr << "#------------------------------------------------------"<< std::endl;
         cfg_file_ofstr << "# Atomistic spin configuration file for vampire"<< std::endl;
         cfg_file_ofstr << "#------------------------------------------------------"<< std::endl;
         cfg_file_ofstr << "# Date: "<< asctime(timeinfo);
         cfg_file_ofstr << "#------------------------------------------------------"<< std::endl;
         cfg_file_ofstr << "Number of spins: "<< vout::total_output_atoms << std::endl;
         cfg_file_ofstr << "System dimensions:" << cs::system_dimensions[0] << "\t" << cs::system_dimensions[1] << "\t" << cs::system_dimensions[2] << std::endl;
         cfg_file_ofstr << "Coordinates-file: atoms-coord.cfg"<< std::endl;
         cfg_file_ofstr << "Time: " << double(sim::time)*mp::dt_SI << std::endl;
         cfg_file_ofstr << "Field: " << sim::H_applied << std::endl;
         cfg_file_ofstr << "Temperature: "<< sim::temperature << std::endl;
         cfg_file_ofstr << "Magnetisation: " << stats::system_magnetization.output_normalized_magnetization() << std::endl;
         cfg_file_ofstr << "Number of Materials: " << mp::num_materials << std::endl;
         for(int mat=0;mat<mp::num_materials;mat++){
            cfg_file_ofstr << mp::material[mat].mu_s_SI << std::endl;
         }
         cfg_file_ofstr << "#------------------------------------------------------" << std::endl;
         cfg_file_ofstr << "Number of spin files: " << vmpi::num_processors-1 << std::endl;
         for(int p=1;p<vmpi::num_processors;p++){
            std::stringstream cfg_sstr;
            cfg_sstr << "atoms-" << std::setfill('0') << std::setw(5) << p << "-" << std::setfill('0') << std::setw(8) << output_atoms_file_counter << ".cfg";
            cfg_file_ofstr << cfg_sstr.str() << std::endl;
         }
         cfg_file_ofstr << "#------------------------------------------------------"<< std::endl;
      }

      // Everyone now outputs their atom list
      cfg_file_ofstr << vout::local_output_atom_list.size() << std::endl;
      for(int i=0; i<vout::local_output_atom_list.size(); i++){
         const int atom = vout::local_output_atom_list[i];
         cfg_file_ofstr << atoms::x_spin_array[atom] << "\t" << atoms::y_spin_array[atom] << "\t" << atoms::z_spin_array[atom] << std::endl;
      }

      cfg_file_ofstr.close();

      output_atoms_file_counter++;

   }