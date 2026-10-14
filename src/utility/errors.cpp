//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2009-2012. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
///
/// @file
/// @brief Contains Vampire error functions and variables
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    10/03/2011
/// @internal
///   Created:    10/03/2011
///   Revision:     ---
///=====================================================================================
///

// Standard Libraries
#include <cstdlib>
#include <iostream>

// vampire header files
#include "errors.hpp"
#include "vmpi.hpp"
#include "vio.hpp"

namespace err
   {
   bool check=false;

   /// Default vexit
   void vexit(){

      // check calling of routine if error checking is activated
      if(err::check==true){
         terminaltextcolor(RED);
         std::cout << "err::vexit has been called" << std::endl;
        terminaltextcolor(WHITE);
      }

      // Abort MPI processes for parallel execution
      #ifdef MPICF
      terminaltextcolor(RED);
      std::cerr << "Fatal error on rank: " << vmpi::my_rank << ": Aborting program." << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << "Fatal error on rank: " << vmpi::my_rank << ": Aborting program." << std::endl;
      // concatenate log and sort
      /*#ifdef WIN_COMPILE
         system("type log.* 2>NUL | sort > log");
      #else
         system("ls log.* | xargs cat | sort -n > log");
      #endif*/

      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      // MPI program dies ungracefully here
      #endif

      // Print error message to screen and log
      zlog << zTs() << "Fatal error: Aborting program." << std::endl;
      terminaltextcolor(RED);
      std::cout << "Fatal error: Aborting program. See log file for details." << std::endl;
      terminaltextcolor(WHITE);

      // Now exit program disgracefully
      exit(EXIT_FAILURE);
   }

   //-----------------------------------------------------------------------------------------
   // Parallel exit function calling MPI_Finalise (must be called by all processes)
   //-----------------------------------------------------------------------------------------
   void v_parallel_all_exit(std::string message){

      #ifdef MPICF

      // Print non-blank error message to screen and log
      std::string blank = "";
      if( message != blank ){
         zlog << zTs() << message << std::endl;
      }
      zlog << zTs() << "Fatal error: Aborting program." << std::endl;
      terminaltextcolor(RED);
      std::cout << "Fatal error: Aborting program. See log file for details." << std::endl;
      terminaltextcolor(WHITE);

      // Wait for everyone
      vmpi::barrier();

      // Now finalise MPI
      vmpi::finalise();
      #endif

      // Now exit program disgracefully
      exit(EXIT_FAILURE);
   }

   //-----------------------------------------------------------------------------------------
   /// zexit with error message
   //-----------------------------------------------------------------------------------------
   void zexit(std::string message){

      // check calling of routine if error checking is activated
      if(err::check==true){std::cout << "err::zexit(std::string) has been called" << std::endl;}

      // Abort MPI processes for parallel execution
      #ifdef MPICF
     terminaltextcolor(RED);
      std::cerr << "Fatal error on rank " << vmpi::my_rank << ": " << message << std::endl;
      std::cerr << "Aborting program." << std::endl;
     terminaltextcolor(WHITE);
      zlog << zTs() << "Fatal error on rank " << vmpi::my_rank << ": " << message << std::endl;
      zlog << zTs() << "Aborting program." << std::endl;

      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      // MPI program dies ungracefully here
      #else
         // Print Error message to screen and log
        terminaltextcolor(RED);
         std::cerr << "Fatal error: " << message << std::endl;
       terminaltextcolor(WHITE);
         zlog << zTs() << "Fatal error: " << message << std::endl;

      #endif

      // Print error message to screen and log
      zlog << zTs() << "Aborting program." << std::endl;
     terminaltextcolor(RED);
      std::cout << "Aborting program. See log file for details." << std::endl;
     terminaltextcolor(WHITE);

      // Now exit program disgracefully
      exit(EXIT_FAILURE);
   }

   //-----------------------------------------------------------------------------------------
   // Fatal error: print message in red to screen and to the log, then exit via vexit.
   // The message is printed as given (callers include "Error - ..." themselves).
   //-----------------------------------------------------------------------------------------
   void fatal(std::string const message){

      terminaltextcolor(RED);
      std::cerr << message << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << message << std::endl;

      err::vexit();

   }

   //-----------------------------------------------------------------------------------------
   // Warning: print message in yellow to screen (rank 0 only) and to the log (every rank),
   // then continue. The message is printed as given (callers include "Warning - ..." themselves).
   //-----------------------------------------------------------------------------------------
   void warning(std::string const message){

      if(vmpi::my_rank == 0){
         terminaltextcolor(YELLOW);
         std::cerr << message << std::endl;
         terminaltextcolor(WHITE);
      }
      zlog << zTs() << message << std::endl;

   }

} // end of namespace err
