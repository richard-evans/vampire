//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Rory Pond 2016. All rights reserved.
//
//   Email: rory.pond@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
// Headers
#include "atoms.hpp"
#include "cells.hpp"
#include "config.hpp"
#include "create.hpp"
#include "demag.hpp"
#include "errors.hpp"
#include "grains.hpp"
#include "ltmp.hpp"
#include "voronoi.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "units.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

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
};