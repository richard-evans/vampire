//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2016, Jack B Collings 2021. All rights reserved.
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
            vin::check_for_valid_value(dl, word, line, prefix, unit, "length", 0.1, 100.0,"input","0.1 - 100.0 Angstroms");
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
         test="RKKYkf";
         if(word==test){
            double kf = atof(value.c_str());
            // Test for valid range
            vin::check_for_valid_value(kf, word, line, prefix, unit, "length", 0.0, 1e19, "input", "0.0 - 100000.0");
            return true;
         }
         //-----------------------
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
            test="material-exponential";
            if(value==test){
               uc::internal::exchange_function = uc::internal::material_exponential;
               return true;
            }
            test="RKKY";
            if(value==test){
               uc::internal::exchange_function = uc::internal::RKKY;
               return true;
            }
            else{
               terminaltextcolor(RED);
               std::cerr << "Error - value for \'exchange:" << word << "\' must be one of:" << std::endl;
               std::cerr << "\t\"nearest-neighbour\"" << std::endl;
               std::cerr << "\t\"shell\"" << std::endl;
               std::cerr << "\t\"exponential\"" << std::endl;
               std::cerr << "\t\"material-exponential\"" << std::endl;
               std::cerr << "\t\"RKKY\"" << std::endl;
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

   // Overloaded match_input_parameter to take in functionalities with super and sub indicies
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line, int superIndex, int subIndex){
      
      std::string prefix = "exchange";
      std::string test = "";

      if(key == prefix){
         test = "unit-cell-category-exchange-parameters";
         if (word == test){
            // Temporary storage container
            std::vector<double> u(3);

            // Read values from string into container
            u = vin::doubles_from_string(value);

            // Check for valid vector
            std::vector <double> range_min(0);  // Holds minimum values for each vector element
            range_min.push_back(0.0);
            range_min.push_back(0.001);
            range_min.push_back(-10000);
            std::vector <double> range_max(0);  // Holds maximum values for each vector element
            range_max.push_back(10000);
            range_max.push_back(100);
            range_max.push_back(10000);

            vin::check_for_valid_vector(u, word, line, prefix, unit, "length", range_min, range_max, "input", "A = 0.0 - 10000, B = 0.001 - 100, C = -10000 - 10000");
            
            // Set values
            internal::material_exchange_parameters[superIndex][subIndex].decay_multiplier = u.at(0);
            internal::material_exchange_parameters[superIndex][subIndex].decay_length = u.at(1);
            internal::material_exchange_parameters[superIndex][subIndex].decay_shift = u.at(2);
            return true;
         }
         test = "nn-cutoff-range";
         if (word == test){
            double nncr=atof(value.c_str());
            vin::check_for_valid_value(nncr, word, line, prefix, unit, "none", 0.0, 1000.0,"input","0.0 - 1000.0 default nearest neighbour distances");
            internal::nn_cutoff_range[superIndex][subIndex] *= nncr;
            return true;
         }
         test = "interaction-cutoff-range";
         if (word == test){
            double icr = atof(value.c_str());
            vin::check_for_valid_value(icr, word, line, prefix, unit, "none", 1.0, 1000.0,"input","1.0 - 1000.0 nearest neighbour distances");
            internal::interaction_cutoff_range[superIndex][subIndex] *= icr;
            return true;
         }
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
