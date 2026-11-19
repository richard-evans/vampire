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
#include <cstdlib>
#include <iostream>
#include <csignal>
#include <sys/wait.h>

// module headers
#include "internal.hpp"

namespace vt{

// helper function overloading standard cstdlib system call to make code cleaner. Every command is stopped after
// TIME_LIMIT seconds (environment variable, default 300, as for the regression tests), so a hang fails its test
// instead of blocking the whole suite.
int system(std::string const &cmd) {

   int limit = 300;
   if(const char* env = std::getenv("TIME_LIMIT")){
      if(std::atoi(env) > 0) limit = std::atoi(env);
   }

   // run the command in its own shell under an alarm; quote it for sh -c
   std::string quoted = "'";
   for(const char c : cmd) quoted += (c == '\'') ? std::string("'\\''") : std::string(1, c);
   quoted += "'";
   const std::string limited = "perl -e 'alarm shift; exec @ARGV' " + std::to_string(limit) + " sh -c " + quoted;

   const int status = std::system(limited.c_str());

   // killed by the alarm: the shells usually hand over to the command, so the status shows SIGALRM itself;
   // a shell that waits for the command instead reports it as exit status 128 + SIGALRM
   const bool timed_out = status != -1 && ( (WIFSIGNALED(status) && WTERMSIG(status) == SIGALRM) ||
                                            (WIFEXITED(status) && WEXITSTATUS(status) == 128 + SIGALRM) );
   if(timed_out){
      std::cerr << "Error: '" << cmd << "' was stopped after " << limit << " s (TIME_LIMIT). Returning as failed test." << std::endl;
   }

   return status;

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

// remove files left by a test run, ignoring any that do not exist
void remove(const std::vector<std::string> files){
   for(const std::string& file : files){
      std::error_code ec;
      std::filesystem::remove(file, ec);
   }
}

}
