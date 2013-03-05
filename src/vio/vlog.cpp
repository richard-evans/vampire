//-----------------------------------------------------------------------------
//
//  Vampire - A code for atomistic simulation of magnetic materials
//
//  Copyright (C) 2009-2012 R.F.L.Evans
//
//  Email:richard.evans@york.ac.uk
//
//  This program is free software; you can redistribute it and/or modify 
//  it under the terms of the GNU General Public License as published by 
//  the Free Software Foundation; either version 2 of the License, or 
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful, but 
//  WITHOUT ANY WARRANTY; without even the implied warranty of 
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
//  General Public License for more details.
//
//  You should have received a copy of the GNU General Public License 
//  along with this program; if not, write to the Free Software Foundation, 
//  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
//
// ----------------------------------------------------------------------------
//
//  Revised log file functions for vampire
//
//  (c) R F L Evans 2013
//
//-----------------------------------------------------------------------------
//
// c++ headers
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// unix headers
#include <time.h>
#include <sys/types.h>
#include <unistd.h>

// vampire headers
#include "errors.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

namespace vio{

// Log output stream
std::ofstream vlog;

// Sundry vlog variables for pretty timestamps
std::string vlog_program_name; // Program Name
std::string vlog_hostname; // Host Name
pid_t       vlog_pid; // Process ID
bool        vlog_is_initialised=false; // Initialised flag

// Function to take command name as argument and initialise timestamp
void initialise_vlog_timestamp(std::string command_name){
   
   // Get program name and process ID
   std::string tmprev;
   int linelength = command_name.length();

   // set character triggers
   const char* key="/"; // Word identifier

   // copy characters after last /
   for(int i=linelength-1;i>=0;i--){

      char c=command_name.at(i);

      if(c != *key){
         tmprev.push_back(c);
      }
      else break;
   }
   
   //reverse read into program name
   linelength=tmprev.size();
   for(int i=linelength-1;i>=0;i--){
      char c=tmprev.at(i);
      vlog_program_name.push_back(c);
   }

   // Get hostname
   char loghostname [80];
   int ghs=gethostname(loghostname, 80);
   if(ghs!=0) std::cerr << "Warning: Unable to retrieve hostname for zlog file." << std::endl; 
   vio::vlog_hostname = loghostname;
   
   // Now get process ID
   vio::vlog_pid = getpid();
         
   // Set unique filename for log if num_procs > 1
   std::stringstream logfn;
   if(vmpi::num_processors==1) logfn << "log2";
   else logfn << "log."<<vmpi::my_rank;
   
   // Open log filename
   std::string log_file = logfn.str();
   const char* log_filec = log_file.c_str();
   vlog.open(log_filec);
   
   // Mark as initialised;
   vlog_is_initialised=true;
   
   vlog << vio::vTs() << "Logfile opened" << std::endl;
   
   return;
}

// Function to generate timestamp for use in vlog fstream
std::string vTs(){

  std::string NullString;
  NullString="";

  if(vio::vlog_is_initialised==true){
     std::ostringstream Ts;
      
     // varibale for time
     time_t seconds;

     // get current time
     seconds = time (NULL);
     struct tm * timeinfo;
     char logtime [80];

     timeinfo = localtime ( &seconds );
     // Format time string
     strftime (logtime,80,"%Y-%m-%d %X ",timeinfo);
  
     Ts << logtime << vio::vlog_program_name << " [" << vio::vlog_hostname << ":" << vio::vlog_pid << ":"<< vmpi::my_rank << "] ";
   
     return Ts.str();
   }
   else{
      std::cerr << "Program error: vlog not initialised. Exiting." << std::endl;
      err::vexit();
   }

   return NullString;
}

} // end of namespace vio

