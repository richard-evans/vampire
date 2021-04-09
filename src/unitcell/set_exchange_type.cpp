//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "exchange.hpp"
#include "errors.hpp"
#include "vio.hpp"

// exchange module headers
#include "internal.hpp"

namespace unitcell{

   //------------------------------------------------------------------------------
   //
   //    Member function of exchange template class
   //
   //    Function to set the type of exchange used in the code from a string
   //    processed from the unit cell file. Default is normalised isotropic
   //    exchange where the values of the exchange are set from the material
   //    file.
   //
   //    Function returns the number of exchange interactions specified in the
   //    unit cell file.
   //
   //------------------------------------------------------------------------------
   unsigned int unitcell::exchange_template_t::set_exchange_type(std::string exchange_type_string){

      //----------------------------------------
      // check for standard isotropic exchange
      //----------------------------------------
      const std::string isotropic_str = "isotropic";
      if(isotropic_str == exchange_type_string){

         // set exchange type
         exchange_type = exchange::isotropic;

         // unset normalization flag
         use_material_exchange_constants = false;

         return 1; // number of exchange interactions

      }

      //----------------------------------------
      // check for standard vectorial exchange
      //----------------------------------------
      const std::string vectorial_str = "vectorial";
      if(vectorial_str == exchange_type_string){

         // set exchange type
         exchange_type = exchange::vectorial;

         // unset normalization flag
         use_material_exchange_constants = false;

         return 3; // number of exchange interactions

      }

      //----------------------------------------
      // check for standard tensorial exchange
      //----------------------------------------
      const std::string tensorial_str = "tensorial";
      if(tensorial_str == exchange_type_string){

         // set exchange type
         exchange_type = exchange::tensorial;

         // unset normalization flag
         use_material_exchange_constants = false;

         return 9; // number of exchange interactions

      }

      //-----------------------------------------
      // check for normalised isotropic exchange
      //-----------------------------------------
      const std::string norm_isotropic_str = "normalised-isotropic";
      if(norm_isotropic_str == exchange_type_string){

         // set exchange type
         exchange_type = exchange::isotropic;

         // unset normalization flag
         use_material_exchange_constants = true;

         return 1; // number of exchange interactions

      }

      //-----------------------------------------
      // check for normalised vectorial exchange
      //-----------------------------------------
      const std::string norm_vectorial_str = "normalised-vectorial";
      if(norm_vectorial_str == exchange_type_string){

         // set exchange type
         exchange_type = exchange::vectorial;

         // unset normalization flag
         use_material_exchange_constants = true;

         return 3; // number of exchange interactions

      }

      //-----------------------------------------
      // check for normalised tensorial exchange
      //-----------------------------------------
      const std::string norm_tensorial_str = "normalised-tensorial";
      if(norm_tensorial_str == exchange_type_string){

         // set exchange type
         exchange_type = exchange::tensorial;

         // unset normalization flag
         use_material_exchange_constants = true;

         return 9; // number of exchange interactions

      }

      terminaltextcolor(RED);
         zlog << zTs() << "\nError: Unknown exchange type \"" << exchange_type_string << "\" in unit cell file. Exiting!" << std::endl;
         std::cerr     << "\nError: Unknown exchange type \"" << exchange_type_string << "\" in unit cell file. Exiting!" << std::endl;
      terminaltextcolor(WHITE);

      // exit program
      err::vexit();

      return 0;

   }

} // end of unitcell namespace
