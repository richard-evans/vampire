//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans and P Chureemart 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <vector>
#include <iomanip>
#include <fstream>
#include <sstream>

// Vampire headers
#include "spintorque.hpp"
#include "vmpi.hpp"
#include "vio.hpp"
#include "sim.hpp"
#include "material.hpp"

// Spin Torque headers
#include "internal.hpp"

namespace st{
   namespace internal{
      //-----------------------------------------------------------------------------
      // Function to output base microcell properties
      //-----------------------------------------------------------------------------
      void output_base_microcell_data(){

         using st::internal::beta_cond;
         using st::internal::beta_diff;
         using st::internal::sa_infinity;
         using st::internal::lambda_sdl;
         using st::internal::pos;

         const int num_cells = beta_cond.size();

         // only output on root process
         if(vmpi::my_rank==0){
            if(sim::time%(ST_output_rate) ==0){
               zlog << zTs() << "Outputting ST base microcell data" << std::endl;
               std::ofstream ofile;
               ofile.open("st-microcells-base.cfg");
               ofile << num_cells << std::endl;
               for(int cell=0; cell < num_cells; ++cell){
	               ofile << pos[3*cell+0] << "\t" << pos[3*cell+1] << "\t" << pos[3*cell+2];
	               ofile << "\t" << beta_cond[cell] << "\t" << beta_diff[cell] << "\t" << sa_infinity[cell] << "\t" << lambda_sdl[cell] << std::endl;
               }

		         ofile.close();
	         }

         }

      return;
      }

      //-----------------------------------------------------------------------------
      // Function to output base microcell properties
      //-----------------------------------------------------------------------------
      void output_microcell_data(){

         using st::internal::m;
         using st::internal::sa;
         using st::internal::j;
         using st::internal::coeff_ast;
         using st::internal::coeff_nast;
         using st::internal::ast;
         using st::internal::nast;
         using st::internal::cell_natom;

         const int num_cells = m.size()/3;

         // only output on root process
         if(vmpi::my_rank==0){

            if(sim::time%(ST_output_rate) ==0){

               // determine file name
               std::stringstream filename;
               filename << "st-microcells-" << std::setfill ('0') << std::setw (8) << config_file_counter << ".cfg";

               zlog << zTs() << "Outputting ST microcell data " << filename.str() << std::endl;

               std::ofstream ofile;
               ofile.open(std::string(filename.str()).c_str());

	            ofile<<"Time:"<< "\t" << sim::time*mp::dt_SI<< std::endl;

               for(int cell=0; cell<num_cells; ++cell){
                  ofile << pos[3*cell+0] << "\t" << pos[3*cell+1] << "\t" << pos[3*cell+2] << "\t";
                  ofile << m[3*cell+0] << "\t" << m[3*cell+1] << "\t" << m[3*cell+2] << "\t";
                  ofile << sa[3*cell+0] << "\t" << sa[3*cell+1] << "\t" << sa[3*cell+2] << "\t";
                  ofile << j[3*cell+0] << "\t" << j[3*cell+1] << "\t" << j[3*cell+2] << "\t";
                  ofile << coeff_ast[cell] << "\t";
                  ofile << coeff_nast[cell] << "\t";
                  ofile << ast[3*cell+0] << "\t" << ast[3*cell+1] << "\t" << ast[3*cell+2] << "\t";
                  ofile << nast[3*cell+0] << "\t" << nast[3*cell+1] << "\t" << nast[3*cell+2] << "\t";
                  ofile << total_ST[3*cell+0] << "\t" << total_ST[3*cell+1] << "\t" << total_ST[3*cell+2];
                  ofile << "\t" << cell_natom[cell];
                  ofile << std::endl;

                  // ofile << spin_torque[3*cell+0] << "\t" << spin_torque[3*cell+1] << "\t" << spin_torque[3*cell+2] << std::endl;
               }

            ofile.close();

            // update config_file_counter
            config_file_counter++;

           }

         }

         return;
      }

   } // end of internal namespace
} // end of st namespace
