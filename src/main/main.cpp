//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2023. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
#include <iostream>
#include <vector>
#include <sstream>

#include "create.hpp"
#include "errors.hpp"
#include "info.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "vmpi.hpp"
#include "vio.hpp"

#include "internal.hpp"

// main namespace
namespace vmain{
   namespace internal{
      std::string input_file_name = "input"; // default input file name
   }
}
int simulate_system();

/// Main function for vampire
/// Prints out program header and calls main program routines
int main(int argc, char* argv[]){
  vout::output_file_name="output"; // default output file name
   // For parallel execution intialise MPI
   vmpi::initialise(argc, argv);

   // Check for valid command-line arguments
   vmain::internal::command_line_args(argc, argv);

   // Initialise log file
   vout::zLogTsInit(std::string(argv[0]));

   // Output Program Header
   if(vmpi::my_rank==0){
      std::cout << "                                                _          " << std::endl;
      std::cout << "                                               (_)         " << std::endl;
      std::cout << "                    __   ____ _ _ __ ___  _ __  _ _ __ ___ " << std::endl;
      std::cout << "                    \\ \\ / / _` | '_ ` _ \\| '_ \\| | '__/ _ \\" << std::endl;
      std::cout << "                     \\ V / (_| | | | | | | |_) | | | |  __/" << std::endl;
      std::cout << "                      \\_/ \\__,_|_| |_| |_| .__/|_|_|  \\___|" << std::endl;
      std::cout << "                                         | |               " << std::endl;
      std::cout << "                                         |_|               " << std::endl;
      std::cout << std::endl;
      std::cout << "                      Version " << vinfo::version() << " " << __DATE__ << " " << __TIME__ << std::endl;
      std::cout << std::endl;
      std::cout << "             Git commit: " << vinfo::githash() << std::endl;
      std::cout << std::endl;
      std::cout << "  Licensed under the GNU Public License(v2). See licence file for details." << std::endl;
      std::cout << std::endl;
      std::cout << "  Developers:   Richard F L Evans, Sarah Jenkins, Andrea Meo, " << std::endl;
      std::cout << "                Daniel Meilak, Andrew Naden, Matthew Ellis," << std::endl;
      std::cout << "                Oscar Arbelaez, Sam Morris, Rory Pond, Weijia Fan," << std::endl;
      std::cout << "                Phanwadee Chureemart, Pawel Sobieszczyk, Joe Barker, " << std::endl;
      std::cout << "                Thomas Ostler, Andreas Biternas, Roy W Chantrell," << std::endl;
      std::cout << "                Wu Hong-Ye, Razvan Ababei, Sam Westmoreland," << std::endl;
      std::cout << "                Milton Persson" << std::endl;
      std::cout << " " << std::endl;
      #ifdef COMP
      std::cout << "                Compiled with:  " << COMP << std::endl;
      #endif
      std::cout << "                Compiler Flags: ";
      #ifdef CUDA
      std::cout << "CUDA ";
      #endif
      #ifdef MPICF
      std::cout << "MPI ";
      #endif
      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << "  Vampire includes a copy of the qhull library from C.B. Barber and The "<< std::endl;
      std::cout << "  Geometry Center and may be obtained via http from www.qhull.org." << std::endl;
      std::cout << std::endl;
      std::cout << "================================================================================" << std::endl;
      time_t rawtime = time(NULL);
      struct tm * timeinfo = localtime(&rawtime);
      std::cout<<asctime(timeinfo);
   }


   #ifdef MPICF
      vmpi::hosts();
   #endif

   #ifdef MPICF
      // nullify non root cout stream
      if(vmpi::my_rank!=0){
         vout::nullify(std::cout);
      }
   #endif

   // Initialise system
   mp::initialise(vmain::internal::input_file_name);

   // Create system
   cs::create();

   // Simulate system
   sim::run();

   // Finalise MPI
   #ifdef MPICF
      vmpi::finalise();
   #endif

   zlog << zTs() << "Simulation ended gracefully." << std::endl;
   terminaltextcolor(GREEN);
   std::cout << "Simulation ended gracefully." << std::endl;
   terminaltextcolor(WHITE);


   return EXIT_SUCCESS;

}
