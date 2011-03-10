#ifndef VIO_H_
#define VIO_H_

#include <fstream>
#include <string>

#include <iostream> 

struct null_streambuf 
: public std::streambuf 
{ 
  void overflow(char c) 
  { 
  } 
}; 

//==========================================================
// Global Output Streams
//==========================================================
extern std::ofstream vinfo;
extern std::ofstream vdp;
extern std::ofstream vmag;

namespace vin{
	extern int read(std::string const);
}

namespace vout{
	
	extern std::vector<unsigned int> file_output_list;
	extern std::vector<unsigned int> screen_output_list;
	
	extern bool output_povray;
	extern int output_povray_rate;

	extern bool output_povray_cells;
	extern int output_povray_cells_rate;
	
	extern void data();

	void redirect(std::ostream& strm, std::string filename);
	void nullify(std::ostream& strm);  

}
#endif /*VIO_H_*/
