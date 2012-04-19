#ifndef VIO_H_
#define VIO_H_

#include <fstream>
#include <string>

#include <iostream> 

#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <vmpi.hpp>

struct null_streambuf 
: public std::streambuf 
{ 
  void overflow(char c) 
  { 
  } 
}; 

// Global Output Streams
extern std::ofstream zinfo;
extern std::ofstream zmag;
extern std::ofstream zgrain;
extern std::ofstream zlog;

// global timestap function
extern std::string zTs();

namespace vin{
	extern int read(std::string const);
}

namespace vout{
	
	extern std::vector<unsigned int> file_output_list;
	extern std::vector<unsigned int> screen_output_list;
	extern std::vector<unsigned int> grain_output_list;
	
	extern int output_grain_rate;
	
	extern bool output_atoms_config;
	extern int output_atoms_config_rate;
	
	extern double atoms_output_min[3];
	extern double atoms_output_max[3];
	
	extern bool output_cells_config;
	extern int output_cells_config_rate;
	
	extern bool output_grains_config;
	extern int output_config_grain_rate;
	
	//extern bool output_povray;
	//extern int output_povray_rate;

	//extern bool output_povray_cells;
	//extern int output_povray_cells_rate;
	
	extern void data();
	extern void config();
	extern void zLogTsInit(std::string);
	
	//extern int pov_file();

	void redirect(std::ostream& strm, std::string filename);
	void nullify(std::ostream& strm);  

}
#endif /*VIO_H_*/
