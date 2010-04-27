#ifndef UNITS_H_
#define UNITS_H_

#include <string>

namespace units{
	extern const double pi;
	
	// conversion functions
	extern int convert(std::string, double&, std::string&);
	//extern int revert(std::string, double, std::string);

}


#endif /* UNITS_H_ */