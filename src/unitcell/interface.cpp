//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <algorithm>
#include <string>

// Vampire headers
#include "unitcell.hpp"
#include "errors.hpp"
#include "vio.hpp"

// unitcell module headers
#include "internal.hpp"

namespace unitcell{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for unitcell module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for first prefix
      std::string prefix="unit-cell";
      std::string test = "";

      // Check for second prefix
      prefix="material";
      if(key == prefix){
         //-------------------------------------------------------------------
   		// Get unit cell filename
   		//-------------------------------------------------------------------
         test="unit-cell-file";
   		if(word==test){
   			std::string ucffile=value;
   			// strip quotes
   			ucffile.erase(std::remove(ucffile.begin(), ucffile.end(), '\"'), ucffile.end());
   			test="";
   			// if filename not blank set ucf file name
   			if(ucffile!=test){
   				//std::cout << matfile << std::endl;
   				uc::internal::unit_cell_filename=ucffile;
   				return true;
   			}
   			else{
   				terminaltextcolor(RED);
   				std::cerr << "Error - empty filename in control statement \'material:" << word << "\' on line " << line << " of input file" << std::endl;
   				terminaltextcolor(WHITE);
   				return false;
   			}
   		}
      }

      // Check for third prefix
      prefix="create";
      if(key == prefix){
         //--------------------------------------------------------------------
         test="crystal-structure";
         if(word==test){
            // Strip quotes
            std::string cs=value;
            cs.erase(std::remove(cs.begin(), cs.end(), '\"'), cs.end());
            uc::internal::crystal_structure=cs;
            return true;
         }
         //--------------------------------------------------------------------
         test="crystal-sublattice-materials";
         if(word==test){
            string t = "true";
            string f = "false";
            if(value==t){
               uc::internal::sublattice_materials = true;
               return true;
            }
            else if(value==f){
               uc::internal::sublattice_materials = false;
               return true;
            }
            else {
               uc::internal::sublattice_materials = true;
               return true;
            }
         }
      }
      // Check for final prefix
      prefix="dimensions";
      if(key == prefix){
         test="unit-cell-size";
         if(word==test){
            double a=atof(value.c_str());
            vin::check_for_valid_value(a, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
            uc::internal::unit_cell_size_x=a;
            uc::internal::unit_cell_size_y=a;
            uc::internal::unit_cell_size_z=a;
            return true;
         }
         //--------------------------------------------------------------------
         test="unit-cell-size-x";
         if(word==test){
            double ax=atof(value.c_str());
            vin::check_for_valid_value(ax, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
            uc::internal::unit_cell_size_x=ax;
            return true;
         }
         //--------------------------------------------------------------------
         test="unit-cell-size-y";
         if(word==test){
            double ay=atof(value.c_str());
            vin::check_for_valid_value(ay, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
            uc::internal::unit_cell_size_y=ay;
            return true;
         }
         //--------------------------------------------------------------------
         test="unit-cell-size-z";
         if(word==test){
            double az=atof(value.c_str());
            vin::check_for_valid_value(az, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
            uc::internal::unit_cell_size_z=az;
            return true;
         }
      }
      // add prefix string
      prefix="exchange";
      if(key == prefix){
         //--------------------------------------------------------------------
         test="interaction-range";
         if(word==test){
            double ir=atof(value.c_str());
            // Test for valid range
            vin::check_for_valid_value(ir, word, line, prefix, unit, "none", 1.0, 1000.0,"input","1.0 - 1000.0 nearest neighbour distances");
            uc::internal::exchange_interaction_range = ir;
            return true;
         }
         //--------------------------------------------------------------------
         test="decay-length";
         if(word==test){
            double dl=atof(value.c_str());
            // Test for valid range
            vin::check_for_valid_value(dl, word, line, prefix, unit, "length", 0.1, 100.0,"input","1.0 - 100.0 Angstroms");
            uc::internal::exchange_decay = dl;
            return true;
         }         //--------------------------------------------------------------------
         test="decay-multiplier";
         if(word==test){
            double dm = atof(value.c_str());
            // Test for valid range
            vin::check_for_valid_value(dm, word, line, prefix, unit, "length", 0.0, 10000.0,"input","0.0 - 10000.0");
            uc::internal::exchange_multiplier = dm;
            return true;
         }
         test="decay-shift";
         if (word==test){
            double ds = atof(value.c_str());
            // Test for valid range
            vin::check_for_valid_value(ds, word, line, prefix, unit, "length", -10000.0, 10000.0,"input","-10000.0 - 10000.0");
            uc::internal::exchange_shift = ds;
            return true;
         }
         test="function";
         if(word==test){
            test="nearest-neighbour";
            if(value==test){
               uc::internal::exchange_function = uc::internal::nearest_neighbour;
               return true;
            }
            test="shell";
            if(value==test){
               uc::internal::exchange_function = uc::internal::shell;
               return true;
            }
            test="exponential";
            if(value==test){
               uc::internal::exchange_function = uc::internal::exponential;
               return true;
            }
            else{
               terminaltextcolor(RED);
               std::cerr << "Error - value for \'exchange:" << word << "\' must be one of:" << std::endl;
               std::cerr << "\t\"nearest-neighbour\"" << std::endl;
               std::cerr << "\t\"exponential\"" << std::endl;
               terminaltextcolor(WHITE);
               err::vexit();
            }
         }  //------------------------------------------------------------
      }
      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index){

      std::string prefix="material:";

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of unitcell namespace
