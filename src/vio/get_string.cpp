//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans and Rory Pond 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <algorithm>
#include <sstream>
#include <string>

// Vampire headers
#include "errors.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "vutil.hpp"

// vio module headers
#include "internal.hpp"

namespace vin {

   //---------------------------------------------------------------------------
   // Function to open file on master process and return string on all
   // processes containing file contents
   //---------------------------------------------------------------------------
   std::string get_string(std::string const filename, std::string source_file_name, int line){

      // boolean variable specifying root process
      bool root = false;
      if(vmpi::my_rank == 0) root = true;

      // number of characters in file (needed by all processors)
      uint64_t length = 0;

      // message buffer to store processed string as characters suitable for MPI Broadcast
      std::vector<char> message(0);

      // Read in file on root
      if (root){

         // Save info to log file
         zlog << zTs() << "   Loading text file " << filename.c_str() << " on root process" << std::endl;

         // instantiate timer to calculate bandwidth
         vutil::vtimer_t timer;

         // ifstream declaration
         std::ifstream inputfile;

         // start the timer
         timer.start();

         // Open file
         inputfile.open(filename.c_str());

         // Check for correct opening
         if(!inputfile.is_open()){
            terminaltextcolor(RED);
            if(line >= 0) std::cerr << "Error opening input file \"" << filename << "\" specified on line " << line << " of " << source_file_name << " : File does not exist or cannot be opened! Exiting!" << std::endl;
            else          std::cerr << "Error opening input file \"" << filename << "\" : File does not exist or cannot be opened! Exiting!" << std::endl;
            terminaltextcolor(WHITE);
            if(line >= 0) zlog << zTs() << "Error opening input file \"" << filename << "\" specified on line " << line << " of " << source_file_name << " : File does not exist or cannot be opened!" << std::endl;
            else          zlog << zTs() << "Error opening input file \"" << filename << "\" : File does not exist or cannot be opened!" << std::endl;
            zlog << zTs() << "If file exists then check file permissions to ensure it is readable by the user." << std::endl;
            err::vexit(); // exit program disgracefully
         }

         // load file directly into std::string
         std::string contents( (std::istreambuf_iterator<char>(inputfile)), std::istreambuf_iterator<char>() );

         // stop the timer
         timer.stop();

         // get total number of characters in file
         length = contents.length();

         // reserve correct amount of storage for message
         message.reserve(length);

         // copy contents to message buffer for broadcast
         std::copy(contents.begin(), contents.end(), std::back_inserter(message));

         // calculate size (MB) and bandwith and save to log file
         const double file_size = length * sizeof(char)*1.0e-6;
         const double bandwidth = file_size / timer.elapsed_time();
         zlog << zTs() << "   File successfully loaded on root process with " << length << " characters and size " << file_size << " [ " << bandwidth << " ] MB/s" << std::endl;

      }

      #ifdef MPICF

         // broadcast string size from root (0) to all processors
         MPI_Bcast(&length, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);

         // save progress to log file
         zlog << zTs() << "   Resizing array on all processes" << std::endl;

         // resize message buffer on all processors other than root
         if(!root) message.resize(length);

         // save progress to log file
         zlog << zTs() << "   Broadcasting file from root process" << std::endl;

         // broadcast message buffer from root (0) to all processors
         vmpi::broadcast(message, 0);
         //MPI_Bcast(&message[0], message.size(), MPI_CHAR, 0, MPI_COMM_WORLD);

      #endif

      // return message array cast to a std::string on all processors
      std::string result_str(message.begin(),message.end());
      return result_str;

   }

} // end of namespace vin
