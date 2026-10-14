//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>
#include <fstream>

// Vampire headers
#include "spintransport.hpp"
#include "vmpi.hpp"

// spintransport module headers
#include "internal.hpp"

namespace spin_transport{

//---------------------------------------------------------------------------------------------------------
// Function to update spin torque field
//---------------------------------------------------------------------------------------------------------
void update(const unsigned int num_local_atoms,            // number of local atoms
            const std::vector<double>& atoms_x_spin_array, // x-spin vector of atoms
            const std::vector<double>& atoms_y_spin_array, // y-spin vector of atoms
            const std::vector<double>& atoms_z_spin_array, // z-spin-vector of atoms
            const std::vector<double>& atoms_m_spin_array  // moment of atoms
   ){

   //-------------------------------------------------------------------------
   // check that module is needed - if not do nothing
   //-------------------------------------------------------------------------
   if( spin_transport::internal::enabled == false ) return;

   //-------------------------------------------------------------------------
   // check that it is time to update
   //-------------------------------------------------------------------------
   if(spin_transport::internal::time_counter != spin_transport::internal::update_rate){
      // if not increment counter and do nothing
      spin_transport::internal::time_counter++;
      return;
   }
   // otherwise reset counter and continue
   else{ spin_transport::internal::time_counter = 1; }

   //---------------------------------------------------------------------------------------------------------
   // update cell magnetizations
   //---------------------------------------------------------------------------------------------------------
   if(spin_transport::internal::sublattice){
      spin_transport::internal::calculate_cell_sublattice_magnetization(num_local_atoms, atoms_x_spin_array, atoms_y_spin_array,
                                                                             atoms_z_spin_array, atoms_m_spin_array);
   }
   else{
      spin_transport::internal::calculate_cell_magnetization(num_local_atoms, atoms_x_spin_array, atoms_y_spin_array,
                                                                  atoms_z_spin_array, atoms_m_spin_array);
   }

   //---------------------------------------------------------------------------------------------------------
   // calculate magnetoresistance
   //---------------------------------------------------------------------------------------------------------
   if(spin_transport::internal::sublattice) spin_transport::internal::calculate_sublattice_resistance();
   else spin_transport::internal::calculate_magnetoresistance();

   //---------------------------------------------------------------------------------------------------------
   // test output of cell-level spin transport data
   //---------------------------------------------------------------------------------------------------------

   /*if( vmpi::my_rank == 0 && !spin_transport::internal::sublattice){
      std::ofstream ofile("stdata.txt");
      for(int i =0; i< spin_transport::internal::total_num_cells; i++){
         const double isat = spin_transport::internal::cell_isaturation[i];
         ofile << spin_transport::internal::cell_position[3*i+0] << "\t" <<
                  spin_transport::internal::cell_position[3*i+1] << "\t" <<
                  spin_transport::internal::cell_position[3*i+2] << "\t" <<
                  spin_transport::internal::cell_magnetization[3*i+0] * isat << "\t" <<
                  spin_transport::internal::cell_magnetization[3*i+1] * isat << "\t" <<
                  spin_transport::internal::cell_magnetization[3*i+2] * isat << "\t" <<
                  spin_transport::internal::cell_spin_torque_fields[3*i+0] << "\t" <<
                  spin_transport::internal::cell_spin_torque_fields[3*i+1] << "\t" <<
                  spin_transport::internal::cell_spin_torque_fields[3*i+2] << "\t" <<
                  spin_transport::internal::cell_resistance[i] << "\t" <<
                  spin_transport::internal::cell_spin_resistance[i] << std::endl;
      }
      ofile.close();
      std::cin.get();
   }

   if( vmpi::my_rank == 0 && spin_transport::internal::sublattice){
      const int num_sublattices = spin_transport::internal::num_sublattices;
      std::ofstream ofile("stdata.txt");
      std::cout << "CRS: " << spin_transport::internal::cell_resistance.size() << "\t" << spin_transport::internal::cell_spin_resistance.size() << std::endl;
      for(uint64_t i =0; i< spin_transport::internal::total_num_cells; i++){
         ofile << spin_transport::internal::cell_position[3*i+0] << "\t" <<
                  spin_transport::internal::cell_position[3*i+1] << "\t" <<
                  spin_transport::internal::cell_position[3*i+2] << "\t"; // <<
                  //spin_transport::internal::cell_resistance[i] << "\t"; // <<
                  //spin_transport::internal::cell_spin_resistance[i] << "\t";

                  for( int sl = 0; sl < num_sublattices; sl++ ){

                     const double isat = spin_transport::internal::cell_sl_isaturation[i][sl];
                     double mix = spin_transport::internal::cell_sl_magnetization_x[i][sl] * isat;
                     double miy = spin_transport::internal::cell_sl_magnetization_y[i][sl] * isat;
                     double miz = spin_transport::internal::cell_sl_magnetization_z[i][sl] * isat;
                     ofile << mix << "\t" << miy << "\t" << miz << "\t" << spin_transport::internal::cell_sl_spin_torque_fields_x[i][sl] << "\t" <<
                                                                           spin_transport::internal::cell_sl_spin_torque_fields_y[i][sl] << "\t" <<
                                                                           spin_transport::internal::cell_sl_spin_torque_fields_z[i][sl] << "\t" <<
                                                                           spin_transport::internal::cell_sl_resistance[i][sl] << "\t" <<
                                                                           spin_transport::internal::cell_sl_spin_resistance[i][sl] << "\t";
                                                                           ;
                  }
                  ofile << std::endl;
      }
      ofile.close();
      std::cin.get();
   }*/

   return;

}

} // end of spin_transport namespace
