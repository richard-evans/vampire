#ifndef VIO_H_
#define VIO_H_

#include <fstream>
#include <string>

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
	
}
		
//==========================================================
// Namespace output
//==========================================================
namespace output{
	extern int output_flags[10];
	extern int output_inc[10];
}

#endif /*VIO_H_*/
