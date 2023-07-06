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
#include "vio.hpp"

// vio module headers
#include "internal.hpp"


void terminaltextcolor(enum textcolor color){
#ifdef WIN_COMPILE
 int fincolor=15;
 if(color==RED) fincolor=12; if(color==GREEN) fincolor=10; if(color==YELLOW) fincolor=14;
 if(color==BLUE) fincolor=9; if(color==PURPLE) fincolor=13;
 SetConsoleTextAttribute(GetStdHandle( STD_OUTPUT_HANDLE ), fincolor);
#else
  std::ostringstream fincolor;
  fincolor<< color;
	std::cout << "\033["<<fincolor.str()<<"m";
	std::cerr << "\033["<<fincolor.str()<<"m";
#endif
}
