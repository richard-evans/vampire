//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) rory.pond 2016. All rights reserved.
//
//   Email: rory.pond@york.ac.uk
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
   int match_input_parameter(std::string const word, std::string const value, std::string const unit, int const line){
   //int match_config          (std::string const word, std::string const value, std::string const unit, int const line){
        std::string prefix="config:";
        //if(key!=prefix) return false;

        // System output config variables
        std::string test="atoms";
        if(word==test){
            internal::output_atoms_config=true;
            return EXIT_SUCCESS;
        }
        //-----------------------------------------
        test="atoms-output-rate";
        if(word==test){
            int i=atoi(value.c_str());
            vin::check_for_valid_int(i, word, line, prefix, 1, 1000000,"input","1 - 1,000,000");
        //  vin::check_for_valid_int(i, word, line, prefix, 1, 1000000,"input","1 - 1,000,000");
            internal::output_atoms_config_rate=i;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="output-new";
        if(word==test){
            internal::output_new = true;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="output-legacy";
        if(word==test){
            internal::output_new = false;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="output-binary";
        if(word==test){
            internal::output_data_format = internal::binary;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="output-text";
        if(word==test){
            internal::output_data_format = internal::text;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="output-mpi-io";
        if(word==test){
            internal::output_data_format = internal::binary;
            internal::output_new = true;
            internal::mpi_io = true;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="output-io-processors";
        if(word==test){
            int x=atoi(value.c_str());
            vin::check_for_valid_int(x, word, line, prefix, 1, vmpi::num_processors,"input","1 - 1,000,000");
            vmpi::num_io_processors=x;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="atoms-minimum-x";
        if(word==test){
            double x=atof(value.c_str());
            vin::check_for_valid_value(x, word, line, prefix, "", "none", 0.0, 1.0,"input","0.0 - 1.0");
            internal::atoms_output_min[0]=x;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="atoms-minimum-y";
        if(word==test){
            double y=atof(value.c_str());
            vin::check_for_valid_value(y, word, line, prefix, "", "none", 0.0, 1.0,"input","0.0 - 1.0");
            internal::atoms_output_min[1]=y;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="atoms-minimum-z";
        if(word==test){
            double z=atof(value.c_str());
            vin::check_for_valid_value(z, word, line, prefix, "", "none", 0.0, 1.0,"input","0.0 - 1.0");
            internal::atoms_output_min[2]=z;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="atoms-maximum-x";
        if(word==test){
            double x=atof(value.c_str());
            vin::check_for_valid_value(x, word, line, prefix, "", "none", 0.0, 1.0,"input","0.0 - 1.0");
            internal::atoms_output_max[0]=x;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="atoms-maximum-y";
        if(word==test){
            double y=atof(value.c_str());
            vin::check_for_valid_value(y, word, line, prefix, "", "none", 0.0, 1.0,"input","0.0 - 1.0");
            internal::atoms_output_max[1]=y;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="atoms-maximum-z";
        if(word==test){
            double z=atof(value.c_str());
            vin::check_for_valid_value(z, word, line, prefix, "", "none", 0.0, 1.0,"input","0.0 - 1.0");
            internal::atoms_output_max[2]=z;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="macro-cells";
        if(word==test){
            internal::output_cells_config=true;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="macro-cells-output-rate";
        if(word==test){
            int i=atoi(value.c_str());
            vin::check_for_valid_int(i, word, line, prefix, 0, 1000000,"input","0 - 1,000,000");
            internal::output_cells_config_rate=i;
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="identify-surface-atoms";
        if(word==test){
            sim::identify_surface_atoms=true;
            return EXIT_SUCCESS;
        }
        //-----------------------------------------
        test="field-range-1-minimum";
        if(word==test){
            double H=atof(value.c_str());
            vin::check_for_valid_value(H, word, line, prefix, unit, "field", -1.e4, 1.0e4,"input","+/- 10,000 T");
            internal::field_output_min_1=H;
            return EXIT_SUCCESS;
        }
        //-----------------------------------------
        test="field-range-1-maximum";
        if(word==test){
            double H=atof(value.c_str());
            vin::check_for_valid_value(H, word, line, prefix, unit, "field", -1.e4, 1.0e4,"input","+/- 10,000 T");
            internal::field_output_max_1=H;
            return EXIT_SUCCESS;
        }
        //-----------------------------------------
        test="field-range-2-minimum";
        if(word==test){
            double H=atof(value.c_str());
            vin::check_for_valid_value(H, word, line, prefix, unit, "field", -1.e4, 1.0e4,"input","+/- 10,000 T");
            internal::field_output_min_2=H;
            return EXIT_SUCCESS;
        }
        //-----------------------------------------
        test="field-range-2-maximum";
        if(word==test){
            double H=atof(value.c_str());
            vin::check_for_valid_value(H, word, line, prefix, unit, "field", -1.e4, 1.0e4,"input","+/- 10,000 T");
            internal::field_output_max_2=H;
            return EXIT_SUCCESS;
        }
        //-----------------------------------------
        else{
        terminaltextcolor(RED);
            std::cerr << "Error - Unknown control statement \'config:"<< word << "\' on line " << line << " of input file" << std::endl;
        terminaltextcolor(WHITE);
            return EXIT_FAILURE;
        }
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
