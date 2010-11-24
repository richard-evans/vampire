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
	
	//extern int scr(std::stringstream);
	extern int pov_file();

  void redirect(std::ostream& strm, std::string filename);
  void nullify(std::ostream& strm);  

}
#endif /*VIO_H_*/
