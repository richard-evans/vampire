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
#include <string>
#include <algorithm>
#include <sstream>

// Vampire headers
#include "errors.hpp"
#include "material.hpp"
#include "vio.hpp"

// material module headers
#include "internal.hpp"

namespace mp{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for material module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      std::string test="";

      //-------------------------------------------------------------------
      // Material file
      //-------------------------------------------------------------------
      std::string prefix="material";
      if(key==prefix){
         test="file";
         if(word==test){
            std::string matfile=value;
            // strip quotes
            matfile.erase(std::remove(matfile.begin(), matfile.end(), '\"'), matfile.end());
            test="";
            if(matfile!=test){
               internal::read_mat_file(matfile,line);
               return true;
            }
            else{
               std::ostringstream message;
               message << "Error - empty filename in control statement \'material:" << word << "\' on line " << line << " of input file";
               err::fatal(message.str());
            }
         }
         return false;
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

      std::string test="";

      //--------------------------------------------------------------------
      test="num-materials";
      if(word==test){
         unsigned int nm = atoi(value.c_str());
         vin::check_for_valid_int(nm, word, line, prefix, 1, 100,"material","1 - 100");
         mp::num_materials=nm;
         return true;
      }
      //--------------------------------------------------------------------
      test="material-name";
      if(word==test){
         mp::read_material[super_index].name=value;
         return true;
      }
      //--------------------------------------------------------------------
      test="damping-constant";
      if(word==test){
         double damping=atof(value.c_str());
         vin::check_for_valid_positive_value(damping, word, line, prefix, unit, "none", 0.0, 10.0,"material","0.0 - 10.0");
         mp::read_material[super_index].alpha=damping;
         return true;
      }
      //--------------------------------------------------------------------
      test="equilibration-damping-constant";
      if(word==test){
         double damping=atof(value.c_str());
         vin::check_for_valid_positive_value(damping, word, line, prefix, unit, "none", 0.0, 10.0,"material","0.0 - 10.0");
         mp::read_material[super_index].alpha_eq=damping;
         return true;
      }
      //--------------------------------------------------------------------
      test="atomic-spin-moment";
      if(word==test){
         double mu_s=atof(value.c_str());
         vin::check_for_valid_positive_value(mu_s, word, line, prefix, unit, "moment", 0.1*9.24e-24, 1e8*9.24e-24,"material","0.1 - 1e8 mu_B");
         mp::read_material[super_index].moment_flag=true;
         mp::read_material[super_index].mu_s_SI=mu_s;
         return true;
      }
      //--------------------------------------------------------------------
      test="relative-gamma";
      if(word==test){
         double gr = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_positive_value(gr, word, line, prefix, unit, "none", 0.01, 100.0,"material"," 0.01 - 100.0");
         mp::read_material[super_index].gamma_rel=gr;
         return true;
      }
      //--------------------------------------------------------------------
      test="material-element";
      if(word==test){
         // Test for 3 characters
         if(value.length()>3){
            std::ostringstream message;
            message << "Error - element identifier on line " << line << " of material file must be a maximum of three characters long";
            err::fatal(message.str());
         }
         else{
            // pad value to be equal to 3 characters
            std::string tmp="   ";
            for(unsigned int i=0;i<3;i++){
               if(i<value.length()){
                  tmp.at(i)=value.at(i);
               }
            }
            mp::read_material[super_index].element=tmp;
            return true;
         }
      }
      //--------------------------------------------------------------------
      test="constrained"; // determines use of alternate integrator
      if(word==test){
         std::string t="true";
         std::string f="false";
         if(value==t){
            mp::read_material[super_index].constrained=true;
            return true;
         }
         else if(value==f){
            mp::read_material[super_index].constrained=false;
            return true;
         }
         else {
            std::ostringstream message;
            message << "Error in input file - material[" << super_index+1 << "]:constrained must be either true or false";
            err::fatal(message.str());
         }
      }
      //--------------------------------------------------------------------
      test="non-magnetic";
      /*
      integer non-magnetic [0]
      The default value is 0 for all materials. Valid values are
      remove, (blank) [same as remove] and keep.
      Value = keep causes the material to be identified as non magnetic
      with all atoms of this type KEPT in the simulation.
      Value = 1 remove/blank causes the material to be identified as non magnetic
      with all atoms of this type REMOVED from the simulation.
      The atomic positions of non-magnetic atoms are saved separately
      with the usual atomic spin configuration for post processing.
      The default value is false for all materials. Valid values are
      true, false or (blank) [same as true].
      */

      if(word==test){
         test="keep";
         // keep all atoms in simulation (for efficient parallelization)
         if(value==test){
            // set flag
            mp::read_material[super_index].non_magnetic = 2;
            return true;
         }
         // delete all atoms in simulation (default)
         test="remove";
         if(value==test){
            // set flag
            mp::read_material[super_index].non_magnetic = 1;
            return true;
         }
         test="";
         if(value==test){
            // set flag
            mp::read_material[super_index].non_magnetic = 1;
            return true;
         }
         else{
            std::ostringstream message;
            message << "Error on line " << line << " of material file - material[" << super_index+1 << "]:" << word << " = " << value << " is not a valid option: remove, keep.";
            err::fatal(message.str());
         }
      }
      //--------------------------------------------------------------------
      test="temperature-rescaling-exponent";
      if(word==test){
         double alpha=atof(value.c_str());
         vin::check_for_valid_value(alpha, word, line, prefix, unit, "none", 0.0, 10.0,"material"," 0.0 - 10.0");
         mp::read_material[super_index].temperature_rescaling_alpha=alpha;
         return true;
      }
      //--------------------------------------------------------------------
      test="temperature-rescaling-curie-temperature";
      if(word==test){
         double Tc=atof(value.c_str());
         vin::check_for_valid_value(Tc, word, line, prefix, unit, "none", 0.0, 10000.0,"material"," 0 - 10000 K");
         mp::read_material[super_index].temperature_rescaling_Tc=Tc;
         return true;
      }
      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of mp namespace
