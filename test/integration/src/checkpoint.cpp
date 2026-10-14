//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// module headers
#include "internal.hpp"

namespace{

// read the last data line of the output file
std::string last_output_line(){
   std::ifstream ifile("output");
   std::string line;
   std::string last_line;
   while(getline(ifile, line)){
      if(!line.empty() && line[0] != '#') last_line = line;
   }
   return last_line;
}

// remove output, log and checkpoint files left by a run
void cleanup(){
   vt::remove({"output", "log"});
   for(const auto& entry : std::filesystem::directory_iterator(".")){
      if(entry.path().extension() == ".chk") vt::remove({entry.path().string()});
   }
}

}

//------------------------------------------------------------------------------
// Test that a simulation saved to a checkpoint and continued from it gives
// exactly the same result as an uninterrupted simulation.
//
//    input-full      uninterrupted simulation
//    input-save      same simulation, stopped part way with sim:save-checkpoint=end
//    input-continue  same simulation, continued with sim:load-checkpoint=continue
//
// The last line of the output file must be identical for both.
//------------------------------------------------------------------------------
bool checkpoint_test(const std::string dir, const std::string executable){

   // get root directory
   std::string path = std::filesystem::current_path();

   // fixed-width output for prettiness
   std::stringstream test_name;
   test_name << "Testing checkpoint restart for " << dir;
   std::cout << std::setw(60) << std::left << test_name.str() << " : " << std::flush;

   // change directory
   if( !vt::chdir(path+"/data/"+dir) ) return false;
   cleanup();

   // run uninterrupted, then saved and continued simulations
   std::string full_line;
   std::string continued_line;
   bool ran = vt::system(executable + " --input-file input-full") == 0;
   if(ran){
      full_line = last_output_line();
      vt::remove({"output", "log"});
      ran = vt::system(executable + " --input-file input-save") == 0 &&
            vt::system(executable + " --input-file input-continue") == 0;
      continued_line = last_output_line();
   }

   // cleanup and return to parent directory
   cleanup();
   if( !vt::chdir(path) ) return false;

   if(!ran){
      std::cout << "FAIL | error running vampire" << std::endl;
      return false;
   }

   if(!full_line.empty() && full_line == continued_line){
      std::cout << "OK" << std::endl;
      return true;
   }
   else{
      std::cout << "FAIL | uninterrupted: " << full_line << "\tcontinued: " << continued_line << std::endl;
      return false;
   }

}
