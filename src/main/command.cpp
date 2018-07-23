//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <string>

// Vampire headers
#include "errors.hpp"
#include "info.hpp"
#include "vmpi.hpp"
#include "vio.hpp"

#include "internal.hpp"

namespace vmain{
namespace internal{

//------------------------------------------------------------
// Function to process command line arguments for vampire
//------------------------------------------------------------
void command_line_args(int argc, char* argv[]){

   // Loop over all arguments
   for(int arg = 1; arg < argc; arg++){

      // convert text to std::string
      std::string sw = argv[arg];

      //-----------------------------
      // version information
      //-----------------------------
      if( sw == "--version" || sw == "-v"){
         if(vmpi::my_rank == 0){
            std::cout << "vampire version " << vinfo::version() << std::endl;
            std::cout << "Githash " << vinfo::githash() << std::endl;
            std::cout << "\nThis is free software; see the source for licence information." << std::endl;
         }
         #ifdef MPICF
            vmpi::barrier();
         	MPI_Finalize();
         #endif
         // exit program
         exit(EXIT_SUCCESS);
         return;
      }

      //-----------------------------
      // input file name
      //-----------------------------
      if(sw=="--input-file"){
         // check number of args not exceeded
         if(arg+1 < argc){
            arg++;
            vmain::internal::input_file_name = string(argv[arg]);
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error - no file specified for \'--input-file\' command line option" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
            return;
         }
      }

      else{
         terminaltextcolor(RED);
         std::cerr << "Error - unknown command line parameter \'" << sw << "\'" << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
         return;
      }
   }

   return;

}

} // end of internal namespace
} // end of main namespace
