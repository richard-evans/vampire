//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2022. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// module headers
#include "internal.hpp"

namespace vt{

// helper function overloading standard cstdlib system call to make code cleaner
int system(std::string const &cmd) {
    return std::system(cmd.c_str());
}

// overloaded chdir function using std::fs and returning true if successful
bool chdir(const std::string path){

   // try changing into path
   try{ std::filesystem::current_path(path); }//setting path
   catch(...){
      std::cerr << "Error changing into directory " << path << ". Returning as failed test." << std::endl;
      return false;
   }

   // otherwise return successfully
   return true;

}

}
