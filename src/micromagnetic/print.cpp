//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <sstream>
#include <fstream>

// Vampire headers
#include "micromagnetic.hpp"
#include "vmpi.hpp"

// micromagnetic module headers
#include "internal.hpp"

namespace micromagnetic{
namespace internal{

//------------------------------------------------------------------------------
// Function to print system parameters for micromagnetic model for all cells
//------------------------------------------------------------------------------
void output_system_parameters(std::vector<double>& pos_and_mom_array,
                              std::vector<double>& x_mag_array,
                              std::vector<double>& y_mag_array,
                              std::vector<double>& z_mag_array){

   std::stringstream fname_ss;
   fname_ss << "parameters-" << vmpi::my_rank << ".txt";
   std::ofstream ofile(fname_ss.str());

   const int num_cells = x_mag_array.size();

   for(int cell=0; cell < num_cells ; cell++){

      const double magn = sqrt(x_mag_array[cell]*x_mag_array[cell] +
                               y_mag_array[cell]*y_mag_array[cell] +
                               z_mag_array[cell]*z_mag_array[cell]);
      ofile << cell << "\t"
            << pos_and_mom_array[4*cell+0] << "\t"
            << pos_and_mom_array[4*cell+1] << "\t"
            << pos_and_mom_array[4*cell+2] << "\t"
            << pos_and_mom_array[4*cell+3] << "\t"
            // normalised magnetization
            << x_mag_array[cell]/pos_and_mom_array[4*cell+3] << "\t"
            << y_mag_array[cell]/pos_and_mom_array[4*cell+3] << "\t"
            << z_mag_array[cell]/pos_and_mom_array[4*cell+3] << "\t"
            << magn/pos_and_mom_array[4*cell+3] << "\t"
            << alpha[cell] << "\t"
            << one_o_chi_perp[cell] << "\t"
            << one_o_chi_para[cell] << "\t"
            << gamma[cell] << "\t"
            << ku[cell] << "\t"
            << ku_x[cell] << "\t"
            << ku_y[cell] << "\t"
            << ku_z[cell] << "\t"
            << ms[cell] << "\t"
            << T[cell] << "\t"
            << Tc[cell] << "\t"
            << m_e[cell] << "\t"
            << alpha_para[cell] << "\t"
            << alpha_perp[cell] << std::endl;
   }

   ofile.close();

}

} // end of internal namespace
} // end of micromagnetic namespace
