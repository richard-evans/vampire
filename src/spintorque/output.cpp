//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <vector>
#include <fstream>

// Vampire headers
#include "spintorque.hpp"
#include "vmpi.hpp"

// Spin Torque headers
#include "internal.hpp"

namespace st{
   namespace internal{
      //-----------------------------------------------------------------------------
      // Function to output microcell properties
      //-----------------------------------------------------------------------------
      void output_microcell_data(){
         
         using st::internal::beta_cond;
         using st::internal::beta_diff;
         using st::internal::sa_infinity;
         using st::internal::lambda_sdl;
         using st::internal::pos;
         
         // only output on root process
         std::cout << vmpi::my_rank << "+++" << std::endl;
         if(vmpi::my_rank==0){
            std::ofstream ofile;
            ofile.open("microcells.cfg");
            
            for(int cell=0; cell<beta_cond.size(); ++cell){
               ofile << pos[3*cell+0] << "\t" << pos[3*cell+1] << "\t" << pos[3*cell+2] << "\t" << beta_cond[cell];
               ofile << "\t" << beta_diff[cell] << "\t" << sa_infinity[cell] << "\t" << lambda_sdl[cell] << std::endl;
            }

            ofile.close();

         }

         return;
      }

   } // end of internal namespace
} // end of st namespace

