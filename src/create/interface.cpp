//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "errors.hpp"
#include "create.hpp"
#include "vio.hpp"

// Internal sim header
#include "internal.hpp"

namespace create{

   //-----------------------------------------------------------------------------
   // Function to process input file parameters for create module
   //-----------------------------------------------------------------------------
   bool match_input_parameter(string const key, string const word, string const value, string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="create";
      if(key!=prefix) return false;

      //----------------------------------
      // Now test for all valid options
      //----------------------------------
      std::string test="faceted-particle";
      if(word==test){
         // check for blank value
         test="";
         if(value == test){
            cs::system_creation_flags[1]=7;
            return true;
         }
         // otherwise require 3 numbers for 100,110 and 111 facet radii
         std::vector<double> u(3);
         u=vin::DoublesFromString(value);
         vin::check_for_valid_three_vector(u, word, line, prefix, "input");
         // check for sensible values
         if(u.at(0) < 1.0 || u.at(0) > 1000.0){
      	   terminaltextcolor(RED);
            std::cerr << "Error: 100 particle facet radius for input variable " << prefix << word << " on line " << line << " of input file must be between 1 and 1000." << std::endl;
            terminaltextcolor(WHITE);
            zlog << zTs() << "Error: 100 particle facet radius for input variable " << prefix << word << " on line " << line << " of input file must be between 1 and 1000." << std::endl;
            err::vexit();
         }
         if(u.at(1) < 1.0 || u.at(1) > 1000.0){
      	   terminaltextcolor(RED);
            std::cerr << "Error: 110 particle facet radius for input variable " << prefix << word << " on line " << line << " of input file must be between 1 and 1000." << std::endl;
            terminaltextcolor(WHITE);
            zlog << zTs() << "Error: 110 particle facet radius for input variable " << prefix << word << " on line " << line << " of input file must be between 1 and 1000." << std::endl;
            err::vexit();
         }
         if(u.at(2) < 1.0 || u.at(2) > 1000.0){
      	   terminaltextcolor(RED);
            std::cerr << "Error: 111 particle facet radius for input variable " << prefix << word << " on line " << line << " of input file must be between 1 and 1000." << std::endl;
            terminaltextcolor(WHITE);
            zlog << zTs() << "Error: 111 particle facet radius for input variable " << prefix << word << " on line " << line << " of input file must be between 1 and 1000." << std::endl;
            err::vexit();
         }
         create::internal::faceted_particle_100_radius = u.at(0);
         create::internal::faceted_particle_110_radius = u.at(1);
         create::internal::faceted_particle_111_radius = u.at(2);
         cs::system_creation_flags[1]=7;
         return true;
      }
      test="cone";
      if(word==test){
         cs::system_creation_flags[1]=8;
         return true;
		}
      // check for truncation factor
      test="cone-angle";
      if(word == test){
         double angle=atof(value.c_str());
         vin::check_for_valid_value(angle, word, line, prefix, unit, "none", 0.1,44.9 ,"input","0.1 - 44.9");
         create::internal::cone_angle=angle;
         return true;
		}

      /*std::string test="slonczewski-spin-polarization-unit-vector";
      if(word==test){
         std::vector<double> u(3);
         u=vin::DoublesFromString(value);
         // Test for valid range
         vin::check_for_valid_unit_vector(u, word, line, prefix, "input");
         // save sanitized unit vector
         sim::internal::slonczewski_spin_polarization_unit_vector = u;
         return true;
      }*/
      //--------------------------------------------------------------------
      // input parameter not found here
      return false;
   }

   //----------------------------------------------------------------------------------
   // material parameter match function
   //----------------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index){

      // add prefix string
      std::string prefix="material:";

      // Check for empty material parameter array and resize
      if(create::internal::mp.size() == 0) create::internal::mp.resize(mp::max_materials);

      //------------------------------------------------------------
      std::string test="alloy-host"; // determines host material
      if(word==test){
         // if this keyword is set, then atoms of this type will be scanned for alloy materials
         create::internal::mp[super_index].alloy_master=true;
         terminaltextcolor(YELLOW);
         std::cout << "Warning: Keyword \'alloy-host\' is deprecated and may be removed in a future release. Please use \'host-alloy\' instead." << std::endl;
         terminaltextcolor(WHITE);
         zlog << zTs() << "Warning: Keyword \'alloy-host\' is deprecated and may be removed in a future release. Please use \'host-alloy\' instead." << std::endl;
         return true;
      }
      //--------------------------------------------------------------------
      test="host-alloy"; // determines host material
      if(word==test){
         // if this keyword is set, then atoms of this type will be scanned for alloy materials
         create::internal::mp[super_index].alloy_master=true;
         // check for type of host alloy
         test=""; // blank (assume homogeneous)
         if(value==test){
            create::internal::mp[super_index].host_alloy_distribution = internal::homogeneous;
            return true;
         }
         test="homogeneous"; // default
         if(value==test){
            create::internal::mp[super_index].host_alloy_distribution = internal::homogeneous;
            return true;
         }
         test="random"; // localised distribution
         if(value==test){
            create::internal::mp[super_index].host_alloy_distribution = internal::random;
            return true;
         }
         test="granular"; // create distribution from intrinsic granular structure
         if(value==test){
            create::internal::mp[super_index].host_alloy_distribution = internal::granular;
            return true;
         }
         // otherwise throw an error
         terminaltextcolor(RED);
         std::cerr << "Error - value for \'material[" << super_index << "]:" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"homogeneous\"" << std::endl;
         std::cerr << "\t\"random\"" << std::endl;
         std::cerr << "\t\"granular\"" << std::endl;
         terminaltextcolor(WHITE);
         zlog << zTs() << "Error - value for \'material[" << super_index << "]:" << word << "\' must be one of:" << std::endl;
         zlog << zTs() << "\t\"homogeneous\"" << std::endl;
         zlog << zTs() << "\t\"random\"" << std::endl;
         zlog << zTs() << "\t\"granular\"" << std::endl;
         err::vexit();

         return true;
      }
      //--------------------------------------------------------------------
      test="host-alloy-smoothness"; // determines host material
      if(word==test){
         // check for smoothness value of host alloy dispersion
         test="standard"; // default
         if(value==test){
            create::internal::mp[super_index].host_alloy_smoothness = 2.0;
            return true;
         }
         test="sharp";
         if(value==test){
            create::internal::mp[super_index].host_alloy_smoothness = 1.0;
            return true;
         }
         test="smooth";
         if(value==test){
            create::internal::mp[super_index].host_alloy_smoothness = 5.0;
            return true;
         }
         // otherwise assume number
         double s=atof(value.c_str());
         vin::check_for_valid_value(s, word, line, prefix, unit, "none", 0.0, 10.0,"material"," 0.0 - 10.0");
         create::internal::mp[super_index].host_alloy_smoothness = s;
         return true;
         // otherwise throw an error (need to check here for conversion of string to number)
         /*else{
            terminaltextcolor(RED);
            std::cerr << "Error - value for \'material[" << super_index << "]:" << word << "\' must be one of:" << std::endl;
            std::cerr << "\t\"standard\"" << std::endl;
            std::cerr << "\t\"sharp\"" << std::endl;
            std::cerr << "\t\"smooth\"" << std::endl;
            std::cerr << "\t<value>" << std::endl;
            zlog << zTs() << "Error - value for \'material[" << super_index << "]:" << word << "\' must be one of:" << std::endl;
            zlog << zTs() << "\t\"standard\"" << std::endl;
            zlog << zTs() << "\t\"sharp\"" << std::endl;
            zlog << zTs() << "\t\"smooth\"" << std::endl;
            zlog << zTs() << "\t\"<value>\"" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }*/
         return true;
      }
      //--------------------------------------------------------------------
      test="host-alloy-scale"; // determines host material
      if(word==test){
         double s=atof(value.c_str());
         vin::check_for_valid_value(s, word, line, prefix, unit, "length", 1, 10000.0,"material"," 0.1 - 1000 nm");
         create::internal::mp[super_index].host_alloy_scale = s;
         return true;
      }
      //--------------------------------------------------------------------
      test="save-host-alloy-distribution"; // saves alloy profile to file
      if(word==test){
         create::internal::mp[super_index].save_host_alloy_profile = true;
         create::internal::mp[super_index].save_file_name = value;
         return true;
      }
      //--------------------------------------------------------------------
      test="alloy-fraction"; // determines %mixing for disordered alloys
      if(word==test){
         double af=atof(value.c_str());
         vin::check_for_valid_value(af, word, line, prefix, unit, "none", 0.0, 1.0,"material"," 0.0 - 1.0");
         create::internal::mp[super_index].slave_material[sub_index].fraction=af;
         return true;
      }
      //--------------------------------------------------------------------
      test="alloy-distribution"; // determines type of alloy distribution in slave
      if(word==test){
         // check for distribution adopted by slave material
         test="native"; // (assumes that of host, default)
         if(value==test){
            create::internal::mp[super_index].slave_material[sub_index].slave_alloy_distribution = internal::native;
            return true;
         }
         test="reciprocal"; // (assumes inverse of that of host)
         if(value==test){
            create::internal::mp[super_index].slave_material[sub_index].slave_alloy_distribution = internal::reciprocal;
            return true;
         }
         test="homogeneous"; // (homogeneous distribution, ignores host distribution)
         if(value==test){
            create::internal::mp[super_index].slave_material[sub_index].slave_alloy_distribution = internal::uniform;
            return true;
         }
         // otherwise throw an error
         terminaltextcolor(RED);
         std::cerr << "Error - value for \'material[" << super_index << "]:" << word << "[" << sub_index << "]\' must be one of:" << std::endl;
         std::cerr << "\t\"native\"" << std::endl;
         std::cerr << "\t\"reciprocal\"" << std::endl;
         std::cerr << "\t\"homogeneous\"" << std::endl;
         terminaltextcolor(WHITE);
         zlog << zTs() << "Error - value for \'material[" << super_index << "]:" << word << "[" << sub_index << "]\' must be one of:" << std::endl;
         zlog << zTs() << "\t\"native\"" << std::endl;
         zlog << zTs() << "\t\"reciprocal\"" << std::endl;
         zlog << zTs() << "\t\"homogeneous\"" << std::endl;
         err::vexit();

         return true;

      }
      //--------------------------------------------------------------------
      test="alloy-variance"; // determines range of alloy fraction in host
      if(word==test){
         // check for type of host alloy
         double v=atof(value.c_str());
         vin::check_for_valid_value(v, word, line, prefix, unit, "none", 0.0, 1.0,"material"," 0.0 - 1.0");
         create::internal::mp[super_index].slave_material[sub_index].variance = v;
         return true;
      }

      //--------------------------------------------------------------------
      // keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of namespace create
