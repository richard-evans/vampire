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
#include <string>

// Vampire headers
#include "config.hpp"
#include "sim.hpp"
#include "errors.hpp"
#include "vio.hpp"

// config module headers
#include "internal.hpp"

namespace config{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for config module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="config";
      if(key!=prefix) return false;

      // System output config variables
      std::string test="atoms";
      if(word==test){
          test="end";
          if(value==test){
              internal::output_atoms_config=true; // Save atomic configuration
              internal::output_atoms_config_continuous=false; // do not save atomic configurations during simulation
              internal::output_atoms_config_end=true; // save atomic configurations at the end of simulation
              return true;
          }
          test="continuous";
          if(value==test){
              internal::output_atoms_config=true; // Save atomic configuration
              internal::output_atoms_config_continuous=true; // save atomic configurations during simulation
              internal::output_atoms_config_end=false; // do not save atomic configurations at the end of simulation
              return true;
          }
          test="";
          if(value==test){
              internal::output_atoms_config=true; // Save atomic configuration
              internal::output_atoms_config_continuous=true; // save atomic configurations during simulation
              internal::output_atoms_config_end=false; // do not save atomic configurations at the end of simulation
              return true;
          }
          else{
              terminaltextcolor(RED);
              std::cerr << "Error - value for \'config:" << word << "\' must be one of:" << std::endl;
              std::cerr << "\t\"\"" << std::endl;
              std::cerr << "\t\"end\"" << std::endl;
              std::cerr << "\t\"continuous\"" << std::endl;
              terminaltextcolor(WHITE);
              err::vexit();
          }
      }
      //-----------------------------------------
      test="atoms-output-rate";
      if(word==test){
         int i=atoi(value.c_str());
         vin::check_for_valid_int(i, word, line, prefix, 1, 1000000,"input","1 - 1,000,000");
      //  vin::check_for_valid_int(i, word, line, prefix, 1, 1000000,"input","1 - 1,000,000");
         internal::output_atoms_config_rate=i;
         return true;
      }
      //
      //-----------------------------------------
      test="sld";
      if(word==test){
         internal::sld_format = true;
         return true;
      }
      //--------------------------------------------------------------------
      test="output-format";
      if(word==test){
         test="binary";
         if(value == test){
            config::internal::format = internal::binary;
            return true;
         }
         test="text";
         if(value == test){
            config::internal::format = internal::text;
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error: Value for \'" << prefix << ":" << word << "\' must be one of:" << std::endl;
            std::cerr << "\t\"binary\"" << std::endl;
            std::cerr << "\t\"text\"" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }
      //--------------------------------------------------------------------
      test="output-mode";
      if(word==test){
         test="legacy";
         if(value == test){
            config::internal::mode = internal::legacy;
            return true;
         }
         test="mpi-io";
         if(value == test){
            config::internal::mode = internal::mpi_io;
            return true;
         }
         test="file-per-process";
         if(value == test){
            config::internal::mode = internal::fpprocess;
            return true;
         }
         test="file-per-node";
         if(value == test){
            config::internal::mode = internal::fpnode;
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error: Value for \'" << prefix << ":" << word << "\' must be one of:" << std::endl;
            std::cerr << "\t\"legacy\"" << std::endl;
            std::cerr << "\t\"mpi-io\"" << std::endl;
            std::cerr << "\t\"file-per-process\"" << std::endl;
            std::cerr << "\t\"file-per-node\"" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }
      //--------------------------------------------------------------------
      test="output-nodes";
      if(word==test){
         int x=atoi(value.c_str());
         vin::check_for_valid_int(x, word, line, prefix, 1, 1000000,"input","1 - 1,000,000");
         if(x > vmpi::num_processors){
            zlog << zTs() << "Warning: Number of parallel output nodes set to " << x << " which is greater than the number of processors (" <<
            vmpi::num_processors << ") used in this simulation. Setting output-nodes to " << vmpi::num_processors << "." << std::endl;
            std::cout << "Warning: Number of parallel output nodes set to " << x << " which is greater than the number of processors (" <<
            vmpi::num_processors << ") used in this simulation. Setting output-nodes to " << vmpi::num_processors << "." << std::endl;
            x = vmpi::num_processors;
         }
         config::internal::num_io_groups = x;
         return true;
      }
      //--------------------------------------------------------------------
      test="atoms-minimum-x";
      if(word==test){
         double x=atof(value.c_str());
         vin::check_for_valid_value(x, word, line, prefix, "", "none", 0.0, 1.0,"input","0.0 - 1.0");
         internal::atoms_output_min[0]=x;
         return true;
      }
      //--------------------------------------------------------------------
      test="atoms-minimum-y";
      if(word==test){
         double y=atof(value.c_str());
         vin::check_for_valid_value(y, word, line, prefix, "", "none", 0.0, 1.0,"input","0.0 - 1.0");
         internal::atoms_output_min[1]=y;
         return true;
      }
      //--------------------------------------------------------------------
      test="atoms-minimum-z";
      if(word==test){
         double z=atof(value.c_str());
         vin::check_for_valid_value(z, word, line, prefix, "", "none", 0.0, 1.0,"input","0.0 - 1.0");
         internal::atoms_output_min[2]=z;
         return true;
      }
      //--------------------------------------------------------------------
      test="atoms-maximum-x";
      if(word==test){
         double x=atof(value.c_str());
         vin::check_for_valid_value(x, word, line, prefix, "", "none", 0.0, 1.0,"input","0.0 - 1.0");
         internal::atoms_output_max[0]=x;
         return true;
      }
      //--------------------------------------------------------------------
      test="atoms-maximum-y";
      if(word==test){
         double y=atof(value.c_str());
         vin::check_for_valid_value(y, word, line, prefix, "", "none", 0.0, 1.0,"input","0.0 - 1.0");
         internal::atoms_output_max[1]=y;
         return true;
      }
      //--------------------------------------------------------------------
      test="atoms-maximum-z";
      if(word==test){
         double z=atof(value.c_str());
         vin::check_for_valid_value(z, word, line, prefix, "", "none", 0.0, 1.0,"input","0.0 - 1.0");
         internal::atoms_output_max[2]=z;
         return true;
      }
      //--------------------------------------------------------------------
      test="macro-cells";
      if(word==test){
          test="end";
          if(value==test){
              internal::output_cells_config=true; // Save atomic configuration
              internal::output_cells_config_continuous=false; // do not save atomic configuration at the end of simulation
              internal::output_cells_config_end=true; // save cells configurations at the end of simulation
              return true;
          }
          test="continuous";
          if(value==test){
              internal::output_cells_config=true; // Save atomic configuration
              internal::output_cells_config_continuous=true; // save atomic configuration during simulation
              internal::output_cells_config_end=false; // do not save cells configurations at the end of simulation
              return true;
          }
          test="";
          if(value==test){
              internal::output_cells_config=true; // Save atomic configuration
              internal::output_cells_config_continuous=true; // save atomic configuration during simulation
              internal::output_cells_config_end=false; // do not save cells configurations at the end of simulation
              return true;
          }
          else{
              terminaltextcolor(RED);
              std::cerr << "Error - value for \'config:" << word << "\' must be one of:" << std::endl;
              std::cerr << "\t\"\"" << std::endl;
              std::cerr << "\t\"end\"" << std::endl;
              std::cerr << "\t\"continuous\"" << std::endl;
              terminaltextcolor(WHITE);
              err::vexit();
          }
      }
      //--------------------------------------------------------------------
      test="macro-cells-output-rate";
      if(word==test){
         int i=atoi(value.c_str());
         vin::check_for_valid_int(i, word, line, prefix, 0, 1000000,"input","0 - 1,000,000");
         internal::output_cells_config_rate=i;
         return true;
      }
      //-------------------------------------------------------------------
      test="identify-surface-atoms";
      if(word==test){
         config::internal::identify_surface_atoms = true;
         return true;
      }
      //-----------------------------------------
      test="field-range-descending-minimum";
      if(word==test){
         double H=atof(value.c_str());
         vin::check_for_valid_value(H, word, line, prefix, unit, "field", -1.e4, 1.0e4,"input","+/- 10,000 T");
         internal::field_output_min_1=H;
         return true;
      }
      //-----------------------------------------
      test="field-range-descending-maximum";
      if(word==test){
         double H=atof(value.c_str());
         vin::check_for_valid_value(H, word, line, prefix, unit, "field", -1.e4, 1.0e4,"input","+/- 10,000 T");
         internal::field_output_max_1=H;
         return true;
      }
      //-----------------------------------------
      test="field-range-ascending-minimum";
      if(word==test){
         double H=atof(value.c_str());
         vin::check_for_valid_value(H, word, line, prefix, unit, "field", -1.e4, 1.0e4,"input","+/- 10,000 T");
         internal::field_output_min_2=H;
         return true;
      }
      //-----------------------------------------
      test="field-range-ascending-maximum";
      if(word==test){
         double H=atof(value.c_str());
         vin::check_for_valid_value(H, word, line, prefix, unit, "field", -1.e4, 1.0e4,"input","+/- 10,000 T");
         internal::field_output_max_2=H;
         return true;
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

      // add prefix string
      std::string prefix="material:";

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of config namespace
