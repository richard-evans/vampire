//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2020. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <fstream>
#include <sstream>
#include <ctime>

// Vampire headers
#include "info.hpp"
#include "vio.hpp"

// vio module headers
#include "internal.hpp"

namespace vout{

void zLogTsInit(std::string tmp){

   // Get program name and process ID
   std::string tmprev;
   int linelength = tmp.length();

   // set character triggers
   const char* key="/";	/// Word identifier

   // copy characters after last /
   for(int i=linelength-1;i>=0;i--){

      char c=tmp.at(i);

      if(c != *key){
         tmprev.push_back(c);
      }
      else break;
   }

   //reverse read into program name
   linelength = tmprev.size();
   for(int i = linelength-1; i>=0; i--){
      char c = tmprev.at(i);
      vout::zLogProgramName.push_back(c);
   }

   // Get hostname
   char loghostname [80];
   #ifdef WIN_COMPILE
      DWORD sizelhn = sizeof ( loghostname );
      int GHS=!GetComputerName(loghostname, &sizelhn); //GetComputerName returns true when retrieves hostname
   #else
      int GHS=gethostname(loghostname, 80);
   #endif
  terminaltextcolor(YELLOW);
   if(GHS!=0) std::cerr << "Warning: Unable to retrieve hostname for zlog file." << std::endl;
  terminaltextcolor(WHITE);
   vout::zLogHostName = loghostname;

   // Now get process ID
   #ifdef WIN_COMPILE
      vout::zLogPid = _getpid();
   #else
      vout::zLogPid = getpid();
   #endif

   // Open log filename
   if(vmpi::my_rank==0) zlog.open("log");

   // Mark as initialised;
   vout::zLogInitialised=true;

   zlog << zTs() << "Logfile opened" << std::endl;

   //------------------------------------
   // Determine current directory
   //------------------------------------
   char directory [256];

   #ifdef WIN_COMPILE
      if(_getcwd(directory, sizeof(directory)) == NULL){
         std::cerr << "Fatal getcwd error in datalog." << std::endl;
     }
   #else
      if(getcwd(directory, sizeof(directory)) == NULL){
         std::cerr << "Fatal getcwd error in datalog." << std::endl;
     }
   #endif

   // write system and version information
   zlog << zTs() << "Executable : " << vout::zLogProgramName << std::endl;
   zlog << zTs() << "Host name  : " << vout::zLogHostName << ":" << std::endl;
   zlog << zTs() << "Directory  : " << directory << std::endl;
   zlog << zTs() << "Process ID : " << vout::zLogPid << std::endl;
   zlog << zTs() << "Version    : " << vinfo::version() << std::endl;
   zlog << zTs() << "Githash    : " << vinfo::githash() << std::endl;

   return;
}

} // end of namespace vout

std::string zTs(){

  std::string NullString;
  NullString="";

	if(vout::zLogInitialised==true){
		std::ostringstream Ts;

		// varibale for time
		time_t seconds;

		// get current time
		seconds = std::time (NULL);
		struct tm * timeinfo;
		char logtime [80];

		timeinfo = localtime ( &seconds );
		// Format time string
		//strftime (logtime,80,"%Y-%m-%d %X ",timeinfo);
      strftime (logtime,80,"%d-%m-%Y [%X] ",timeinfo);

		Ts << logtime;
      // << vout::zLogProgramName << " [" << vout::zLogHostName << ":" << vout::zLogPid << ":"<< vmpi::my_rank << "] ";

		return Ts.str();

	}
	else{
		terminaltextcolor(RED);
		std::cerr << "Error! - zlog not initialised, exiting" << std::endl;
		// This can be recursive - vexit calls zTs()
		//err::vexit();
		// Exit manually
		std::cerr << "Fatal error: Aborting program. See log file for details." << std::endl;
	  terminaltextcolor(WHITE);
		exit(EXIT_FAILURE);
	}

	return NullString;
}
