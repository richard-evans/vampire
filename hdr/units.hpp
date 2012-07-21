#ifndef UNITS_H_
#define UNITS_H_

#include <string>
#include <vector>

namespace units{
	extern const double pi;
	
	// conversion functions
	extern int convert(std::string, double&, std::string&);
	extern void convert(std::string, std::vector<double>&, std::string&);

}


#endif /* UNITS_H_ */

