//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>

// Vampire headers
#include "anisotropy.hpp"
#include "cells.hpp"
#include "config.hpp"
#include "create.hpp"
#include "dipole.hpp"
#include "environment.hpp"
#include "exchange.hpp"
#include "grains.hpp"
#include "hierarchical.hpp"
#include "integrate.hpp"
#include "material.hpp"
#include "micromagnetic.hpp"
#include "montecarlo.hpp"
#include "program.hpp"
#include "sim.hpp"
#include "sld.hpp"
#include "spininitialize.hpp"
#include "spintorque.hpp"
#include "spintransport.hpp"
#include "spinwaves.hpp"
#include "stats.hpp"
#include "unitcell.hpp"
#include "vio.hpp"

// material module headers
#include "internal.hpp"

namespace mp{

namespace internal{

   //---------------------------------------------------------------------------
   // Function to match material file keywords to module material parameters
   //---------------------------------------------------------------------------
   int match_material(std::string const word,
                      std::string const value,
                      std::string const unit,
                      int const line,
                      int const super_index,
                      int const sub_index,
                      std::string const line_string,
                      std::string const filename_string)
   {

      //-------------------------------------------------------------------
      // Call module material parameters
      //-------------------------------------------------------------------
      if(mp::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(anisotropy::match_material_parameter(word, value, unit, line, super_index, sub_index, mp::max_materials)) return EXIT_SUCCESS;
      else if(cells::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(config::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(create::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(grains::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(dipole::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(exchange::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(hierarchical::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(integrate::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(montecarlo::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(program::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(sim::match_material_parameter(word, value, unit, line, super_index)) return EXIT_SUCCESS;
      else if(stats::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(st::match_material(word, value, unit, line, super_index)) return EXIT_SUCCESS;
      else if(spin_transport::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(unitcell::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(micromagnetic::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(environment::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(sld::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(spininitialize::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(spinwaves::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;

      //--------------------------------------------------------------------
      // keyword not found
      //--------------------------------------------------------------------
      terminaltextcolor(RED);
      std::cerr << "Error - Unknown control statement '" << line_string << "' on line " << line << " of material file '" << filename_string << "'" << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << "Error - Unknown control statement '" << line_string << " on line " << line << " of material file '" << filename_string << "'" << std::endl;
      return EXIT_FAILURE;

   }


} // end of internal namespace

} // end of mp namespace
