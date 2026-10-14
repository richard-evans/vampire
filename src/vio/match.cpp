//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans and Rory Pond 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//          rory.pond@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <algorithm>
#include <sstream>
#include <iostream>

// Vampire headers
#include "anisotropy.hpp"
#include "cells.hpp"
#include "config.hpp"
#include "create.hpp"
#include "dipole.hpp"
#include "environment.hpp"
#include "errors.hpp"
#include "exchange.hpp"
#include "gpu.hpp"
#include "grains.hpp"
#include "hamr.hpp"
#include "hierarchical.hpp"
#include "integrate.hpp"
#include "ltmp.hpp"
#include "material.hpp"
#include "micromagnetic.hpp"
#include "montecarlo.hpp"
#include "program.hpp"
#include "sim.hpp"
#include "sld.hpp"
#include "spininitialize.hpp"
#include "spintorque.hpp"
#include "spintransport.hpp"
#include "spinwaves.hpp" // JRH
#include "stats.hpp"
#include "unitcell.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// vio module headers
#include "internal.hpp"

namespace vin{

   //---------------------------------------------------------------------------
   // Function to match input file keywords to module input parameters. Each
   // module processes its own keywords in src/<module>/interface.cpp
   //---------------------------------------------------------------------------
   int match(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      //-------------------------------------------------------------------
      // Call module input parameters
      //-------------------------------------------------------------------
      if(ltmp::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(anisotropy::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(cells::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(config::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(create::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(grains::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(dipole::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(gpu::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(exchange::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(hierarchical::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(montecarlo::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(program::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(integrate::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(sim::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(stats::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(st::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(spin_transport::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(unitcell::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(mp::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(vio::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(vmpi::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(micromagnetic::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(environment::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(hamr::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(sld::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(spininitialize::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(spinwaves::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS; // JRH spinwaves input parameters

      //--------------------------------------------------------------------
      // keyword not found
      //--------------------------------------------------------------------
      terminaltextcolor(RED);
      std::cerr << "Error - Unknown control statement \'" << key <<":"<< word << "\' on line " << line << " of input file" << std::endl;
      terminaltextcolor(WHITE);
      return EXIT_FAILURE;

   } // end of match function

   int match(std::string const key, std::string const word, std::string const value, std::string const unit, int const line, int const superIndex, int const subIndex){
      if(unitcell::match_input_parameter(key, word, value, unit, line, superIndex, subIndex)) return EXIT_SUCCESS;

      terminaltextcolor(RED);
      std::cerr << "Error - Unknown control statement \'" << key <<":"<< word << "\' on line " << line << " of input file" << std::endl;
      terminaltextcolor(WHITE);
      return EXIT_FAILURE;
   }
}
