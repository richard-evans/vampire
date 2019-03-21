//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans and Rory Pond 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk and rory.pond@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <sstream>

// Vampire headers
#include "config.hpp"
#include "errors.hpp"
#include "info.hpp"
#include "gpu.hpp"
#include "grains.hpp"
#include "sim.hpp"
#include "vio.hpp"

// vio module headers
#include "internal.hpp"

namespace vout{
	std::string zLogProgramName; /// Program Name
	std::string zLogHostName; /// Host Name
	bool        zLogInitialised=false; /// Initialised flag
	#ifdef WIN_COMPILE
		int      zLogPid; /// Process ID
	#else
		pid_t    zLogPid; /// Process ID
	#endif
}

///-------------------------------------------------------
/// Function to write header information about simulation
///-------------------------------------------------------
void write_output_file_header(std::ofstream& ofile, std::vector<unsigned int>& file_output_list){

	//------------------------------------
	// Determine current time
	//------------------------------------
	time_t seconds;

	// get time now
	seconds = time (NULL);
	struct tm * timeinfo;
	char oftime [80];

	timeinfo = localtime ( &seconds );
	// format time string
	strftime (oftime,80,"%Y-%m-%d %X ",timeinfo);

	//------------------------------------
	// Determine current directory
	//------------------------------------
	char directory [256];

	#ifdef WIN_COMPILE
		_getcwd(directory, sizeof(directory));
	#else
		getcwd(directory, sizeof(directory));
	#endif

	//------------------------------------
	// Output output file header
	//------------------------------------
	ofile << "#----------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
	ofile << "# " << "Output file for vampire simulation" << std::endl;
	ofile << "# " << "  time       : " << oftime << "    process id : " << vout::zLogPid << std::endl;
	ofile << "# " << "  hostname   : " << vout::zLogHostName << std::endl;
	ofile << "# " << "  path       : " << directory << std::endl;
   ofile << "# " << "  version    : " << vinfo::version() << std::endl;
   ofile << "# " << "  githash    : " << vinfo::githash() << std::endl;
	ofile << "#----------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
	//ofile << "# time" << "\t" << "temperature" << "\t" <<  "|m|" << "\t" << "..." << std::endl; // to be concluded...

	return;

}

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
		linelength=tmprev.size();
		for(int i=linelength-1;i>=0;i--){
			char c=tmprev.at(i);
			zLogProgramName.push_back(c);
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
		zLogHostName = loghostname;

		// Now get process ID
		#ifdef WIN_COMPILE
			zLogPid = _getpid();
		#else
			zLogPid = getpid();
		#endif

		// Open log filename
		if(vmpi::my_rank==0) zlog.open("log");

		// Mark as initialised;
		zLogInitialised=true;

		zlog << zTs() << "Logfile opened" << std::endl;

      //------------------------------------
   	// Determine current directory
   	//------------------------------------
   	char directory [256];

   	#ifdef WIN_COMPILE
   		_getcwd(directory, sizeof(directory));
   	#else
   		getcwd(directory, sizeof(directory));
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

	// Data output wrapper function
	void data(){

		// check calling of routine if error checking is activated
		if(err::check==true){std::cout << "vout::data has been called" << std::endl;}

		// Calculate MPI Timings since last data output
		#ifdef MPICF
		if(vmpi::DetailedMPITiming){

			// Calculate Average times
			MPI_Reduce (&vmpi::TotalComputeTime,&vmpi::AverageComputeTime,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			MPI_Reduce (&vmpi::TotalWaitTime,&vmpi::AverageWaitTime,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			vmpi::AverageComputeTime/=double(vmpi::num_processors);
			vmpi::AverageWaitTime/=double(vmpi::num_processors);

			// Calculate Maximum times
			MPI_Reduce (&vmpi::TotalComputeTime,&vmpi::MaximumComputeTime,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
			MPI_Reduce (&vmpi::TotalWaitTime,&vmpi::MaximumWaitTime,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

			// Save times for timing matrix
			vmpi::ComputeTimeArray.push_back(vmpi::TotalComputeTime);
			vmpi::WaitTimeArray.push_back(vmpi::TotalWaitTime);

			// reset until next data output
			vmpi::TotalComputeTime=0.0;
			vmpi::TotalWaitTime=0.0;
		}
		#endif

      // check for open ofstream on root process only
      if(vmpi::my_rank == 0){
         if(!zmag.is_open()){
            // check for checkpoint continue and append data
            if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag) zmag.open("output",std::ofstream::app);
            // otherwise overwrite file
            else{
               zmag.open("output",std::ofstream::trunc);
               // write file header information
               write_output_file_header(zmag, file_output_list);
            }
         }
      }

		// Only output 1/output_rate time steps
		if(sim::time%vout::output_rate==0){

		// Output data to output
		if(vmpi::my_rank==0){

         // For gpu acceleration get statistics from device
         if(gpu::acceleration) gpu::stats::get();

			for(unsigned int item=0;item<file_output_list.size();item++){
				switch(file_output_list[item]){
					case 0:
						vout::time(zmag);
						break;
					case 1:
						vout::real_time(zmag);
						break;
					case 2:
						vout::temperature(zmag);
						break;
					case 3:
						vout::Happ(zmag);
						break;
					case 4:
						vout::Hvec(zmag);
						break;
					case 5:
						vout::mvec(zmag);
						break;
					case 6:
						vout::magm(zmag);
						break;
					case 7:
						vout::mean_magm(zmag);
						break;
					case 8:
						vout::mat_mvec(zmag);
						break;
					case 9:
						vout::mat_mean_magm(zmag);
						break;
					case 12:
						vout::mdoth(zmag);
						break;
					case 14:
						vout::systorque(zmag);
						break;
					case 15:
						vout::mean_systorque(zmag);
						break;
					case 16:
						vout::constraint_phi(zmag);
						break;
					case 17:
						vout::constraint_theta(zmag);
						break;
					case 18:
						vout::material_constraint_phi(zmag);
						break;
					case 19:
						vout::material_constraint_theta(zmag);
						break;
					case 20:
						vout::material_mean_systorque(zmag);
						break;
					case 21:
						vout::mean_system_susceptibility(zmag);
						break;
					case 22:
						vout::phonon_temperature(zmag);
						break;
					case 23:
						vout::material_temperature(zmag);
						break;
					case 24:
						vout::material_applied_field_strength(zmag);
						break;
					case 25:
						vout::material_fmr_field_strength(zmag);
						break;
					case 26:
						vout::mat_mdoth(zmag);
						break;
					case 27:
						vout::total_energy(zmag);
						break;
					case 28:
						vout::mean_total_energy(zmag);
						break;
					case 29:
						vout::total_anisotropy_energy(zmag);
						break;
					case 30:
						vout::mean_total_anisotropy_energy(zmag);
						break;
					case 31:
						//vout::total_cubic_anisotropy_energy(zmag);
						break;
					case 32:
						//vout::mean_total_cubic_anisotropy_energy(zmag);
						break;
					case 33:
						//vout::total_surface_anisotropy_energy(zmag);
						break;
					case 34:
						//vout::mean_total_surface_anisotropy_energy(zmag);
						break;
					case 35:
						vout::total_exchange_energy(zmag);
						break;
					case 36:
						vout::mean_total_exchange_energy(zmag);
						break;
					case 37:
						vout::total_applied_field_energy(zmag);
						break;
					case 38:
						vout::mean_total_applied_field_energy(zmag);
						break;
					case 39:
						vout::total_magnetostatic_energy(zmag);
						break;
					case 40:
						vout::mean_total_magnetostatic_energy(zmag);
						break;
					case 41:
						//vout::total_so_anisotropy_energy(zmag);
						break;
					case 42:
						//vout::mean_total_so_anisotropy_energy(zmag);
						break;
					case 43:
						vout::height_mvec(zmag);
						break;
					case 44:
						vout::material_height_mvec(zmag);
						break;
					case 45:
						vout::height_mvec_actual(zmag);
						break;
					case 46:
						vout::material_height_mvec_actual(zmag);
						break;
					case 47:
						vout::fmr_field_strength(zmag);
						break;
               case 48:
						vout::mean_mvec(zmag);
						break;
               case 49:
						vout::mat_mean_mvec(zmag);
						break;
               case 50:
						vout::mean_material_susceptibility(zmag);
						break;
					case 51:
						vout::mean_height_magnetisation_length(zmag);
						break;
					case 52:
						vout::mean_height_magnetisation(zmag);
						break;
					case 60:
						vout::MPITimings(zmag);
						break;
					case 61:
						vout::mean_system_specific_heat(zmag);
						break;
					case 62:
						vout::mean_material_specific_heat(zmag);
						break;
               case 63:
                  vout::material_total_energy(zmag);
                  break;
               case 64:
                  vout::material_mean_total_energy(zmag);
                  break;
               case 65:
                  vout::resistance(zmag);
                  break;
               case 66:
                  vout::current(zmag);
                  break;
               default:
                  std::cerr     << "\nProgrammer error: unknown output function " << file_output_list[item] << ". Exiting" << std::endl;
                  zlog << zTs() << "Programmer error: unknown output function " << file_output_list[item] << ". Exiting" << std::endl;
                  err::vexit();
				}
			}
			// Carriage return
			if(file_output_list.size()>0) zmag << std::endl;

			} // end of code for rank 0 only
		} // end of if statement for output rate

		// Output data to cout
		if(vmpi::my_rank==0){
			if(sim::time%vout::output_rate==0){ // needs to be altered to separate variable at some point

            // For gpu acceleration get statistics from device (repeated here so additional performance cost for doing this)
            if(gpu::acceleration) gpu::stats::get();

				for(unsigned int item=0;item<screen_output_list.size();item++){
				switch(screen_output_list[item]){
					case 0:
						vout::time(std::cout);
						break;
					case 1:
						vout::real_time(std::cout);
						break;
					case 2:
						vout::temperature(std::cout);
						break;
					case 3:
						vout::Happ(std::cout);
						break;
					case 4:
						vout::Hvec(std::cout);
						break;
					case 5:
						vout::mvec(std::cout);
						break;
					case 6:
						vout::magm(std::cout);
						break;
					case 7:
						vout::mean_magm(std::cout);
						break;
					case 8:
						vout::mat_mvec(std::cout);
						break;
					case 9:
						vout::mat_mean_magm(std::cout);
						break;
					case 12:
						vout::mdoth(std::cout);
						break;
					case 14:
						vout::systorque(std::cout);
						break;
					case 15:
						vout::mean_systorque(std::cout);
						break;
					case 16:
						vout::constraint_phi(std::cout);
						break;
					case 17:
						vout::constraint_theta(std::cout);
						break;
					case 18:
						vout::material_constraint_phi(std::cout);
						break;
					case 19:
						vout::material_constraint_theta(std::cout);
						break;
					case 20:
						vout::material_mean_systorque(std::cout);
						break;
					case 21:
						vout::mean_system_susceptibility(std::cout);
						break;
					case 22:
						vout::phonon_temperature(std::cout);
						break;
					case 23:
						vout::material_temperature(std::cout);
						break;
					case 24:
						vout::material_applied_field_strength(std::cout);
						break;
					case 25:
						vout::material_fmr_field_strength(std::cout);
						break;
					case 26:
						vout::mat_mdoth(std::cout);
						break;
					case 27:
						vout::total_energy(std::cout);
						break;
					case 28:
						vout::mean_total_energy(std::cout);
						break;
					case 29:
						vout::total_anisotropy_energy(std::cout);
						break;
					case 30:
						vout::mean_total_anisotropy_energy(std::cout);
						break;
					case 31:
						//vout::total_cubic_anisotropy_energy(std::cout);
						break;
					case 32:
						//vout::mean_total_cubic_anisotropy_energy(std::cout);
						break;
					case 33:
						//vout::total_surface_anisotropy_energy(std::cout);
						break;
					case 34:
						//vout::mean_total_surface_anisotropy_energy(std::cout);
						break;
					case 35:
						vout::total_exchange_energy(std::cout);
						break;
					case 36:
						vout::mean_total_exchange_energy(std::cout);
						break;
					case 37:
						vout::total_applied_field_energy(std::cout);
						break;
					case 38:
						vout::mean_total_applied_field_energy(std::cout);
						break;
					case 39:
						vout::total_magnetostatic_energy(std::cout);
						break;
					case 40:
						vout::mean_total_magnetostatic_energy(std::cout);
						break;
					case 41:
						//vout::total_so_anisotropy_energy(std::cout);
						break;
					case 42:
						//vout::mean_total_so_anisotropy_energy(std::cout);
						break;
					case 47:
						vout::fmr_field_strength(std::cout);
						break;
               case 48:
						vout::mean_mvec(std::cout);
						break;
               case 49:
						vout::mat_mean_mvec(std::cout);
						break;
               case 50:
						vout::mean_material_susceptibility(std::cout);
						break;
					case 60:
						vout::MPITimings(std::cout);
						break;
					case 61:
						vout::mean_system_specific_heat(std::cout);
						break;
					case 62:
						vout::mean_material_specific_heat(std::cout);
						break;
               case 63:
                  vout::material_total_energy(std::cout);
                  break;
               case 64:
                  vout::material_mean_total_energy(std::cout);
                  break;
               case 65:
                  vout::resistance(std::cout);
                  break;
               case 66:
                  vout::current(std::cout);
                  break;
               default:
                  std::cerr     << "\nProgrammer error: unknown screen output function " << file_output_list[item] << ". Exiting" << std::endl;
                  zlog << zTs() << "Programmer error: unknown screen output function " << file_output_list[item] << ". Exiting" << std::endl;
                  err::vexit();

				}
			}

			// Carriage return
			if(screen_output_list.size()>0) std::cout << std::endl;
			}

		} // End of if statement to output data to screen

		if(sim::time%vout::output_grain_rate==0){

		// calculate grain magnetisations
		grains::mag();

		// Output data to zgrain
		if(vmpi::my_rank==0){

			// check for open ofstream
			if(vout::grain_output_list.size() > 0 && !zgrain.is_open()){
				// check for checkpoint continue and append data
				if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag) zgrain.open("grain",std::ofstream::app);
				// otherwise overwrite file
				else zgrain.open("grain",std::ofstream::trunc);
			}

			for(unsigned int item=0;item<vout::grain_output_list.size();item++){
			switch(vout::grain_output_list[item]){
				case 0:
					vout::time(zgrain);
					break;
				case 1:
					vout::real_time(zgrain);
					break;
				case 2:
					vout::temperature(zgrain);
					break;
				case 3:
					vout::Happ(zgrain);
					break;
				case 4:
					vout::Hvec(zgrain);
					break;
				case 10:
					vout::grain_mvec(zgrain);
					break;
				case 11:
					vout::grain_magm(zgrain);
					break;
				case 13:
					vout::grain_mat_mvec(zgrain);
					break;
				case 22:
					vout::phonon_temperature(zgrain);
					break;
			}
		}

		// Carriage return
		if(vout::grain_output_list.size()>0) zgrain << std::endl;
		}
		}

		// Output configuration files to disk
		config::output();

		// optionally save checkpoint file
		if(sim::save_checkpoint_flag==true && sim::save_checkpoint_continuous_flag==true && sim::time%sim::save_checkpoint_rate==0) save_checkpoint();

	} // end of data
}

/// @brief Function to output timestamp to stream
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2012. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    19/04/2012
///
/// @return TS
///
/// @internal
///	Created:		19/04/2012
///	Revision:	  ---
///=====================================================================================
///
std::string zTs(){

  std::string NullString;
  NullString="";

	if(vout::zLogInitialised==true){
		std::ostringstream Ts;

		// varibale for time
		time_t seconds;

		// get current time
		seconds = time (NULL);
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
