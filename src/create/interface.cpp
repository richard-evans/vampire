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
      test="voronoi-grain-substructure";
      if(word==test){
         create::internal::generate_voronoi_substructure = true;
         return true;
      }
      test="voronoi-grain-substructure-crystallization-radius";
      if(word==test){
         double rsize=atof(value.c_str());
         vin::check_for_valid_value(rsize, word, line, prefix, unit, "none", 0.01, 2.0,"input","0.01 - 2");
         create::internal::voronoi_grain_substructure_crystallization_radius=rsize;
         return true;
      }
      test="voronoi-grain-size";
      if(word==test){
         double psize=atof(value.c_str());
         vin::check_for_valid_value(psize, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
         create::internal::voronoi_grain_size=psize;
         return true;
      }
      else
      //--------------------------------------------------------------------
      test="voronoi-grain-spacing";
      if(word==test){
         double pspacing=atof(value.c_str());
         vin::check_for_valid_value(pspacing, word, line, prefix, unit, "length", 0.0, 1.0e7,"input","0.0 Angstroms - 1 millimetre");
         create::internal::voronoi_grain_spacing=pspacing;
         return true;
      }
      //--------------------------------------------------------------------
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
      //--------------------------------------------------------------------
      test="bubble";
      if(word==test){
         cs::system_creation_flags[1]=9;
         return true;
		}
      //--------------------------------------------------------------------
      test="bubble-radius";
      if(word==test){
         double r=atof(value.c_str());
         vin::check_for_valid_value(r, word, line, prefix, unit, "none", 0.0,1.0 ,"input","0.0 - 1.0");
         create::internal::bubble_radius=r;
         return true;
		}
      //--------------------------------------------------------------------
      test="bubble-nucleation-height";
      if(word==test){
         double nh=atof(value.c_str());
         vin::check_for_valid_value(nh, word, line, prefix, unit, "none", 0.0,1.0 ,"input","0.0 - 1.0");
         create::internal::bubble_nucleation_height=nh;
         return true;
      }
      //--------------------------------------------------------------------
      test="select-material-by-height";
      if(word==test){
          create::internal::select_material_by_z_height = true; // default
          // also check for value
          std::string VFalse="false";
          if(value==VFalse){
             create::internal::select_material_by_z_height = false; // default
          }
          return EXIT_SUCCESS;
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
      if(create::internal::mp.size() == 0){
         create::internal::mp.resize(mp::max_materials);
         // initialise unit cell/material associations. Value should be zero for unit cells with 1 material so that by default CSG operations are applied to those atoms
         for(int i = 0; i < mp::max_materials; i++) create::internal::mp[i].unit_cell_category = 0;
      }

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
      test="fill-substructure-space";
      if(word==test){
         // Test for sane input
         bool sanitised_bool = vin::check_for_valid_bool(value, word, line, prefix,"material");
         // set flag
         create::internal::mp[super_index].sub_fill = sanitised_bool;
         return true;
      }
      /*
         Float to set the reduced starting height (as a fraction of the total system height) for
         the voronoi grain substructure. At this height the voronoi grain size is the standard size.
         Away from the nulceation height the voronoi grain size is reduced according to
         size = (1-x/max)**radius.
      */
      test="voronoi-grain-substructure-nucleation-height";
      if(word==test){
         double nh=atof(value.c_str());
         vin::check_for_valid_value(nh, word, line, prefix, unit, "none", 0.0, 1.0,"material"," 0.0 - 1.0");
         create::internal::mp[super_index].voronoi_grain_substructure_nucleation_height = nh;
         return true;
      }
      /*
         integer to associate the material to a particular material within the unit cell.
         Default is 0 but can be overidden with this parameter.
      */
      test="unit-cell-category";
      if(word==test){
         int uccat=atoi(value.c_str());
         vin::check_for_valid_int(uccat, word, line, prefix, 0, mp::max_materials,"material"," 1 - 100");
         create::internal::mp[super_index].unit_cell_category = uccat - 1; // subtract 1 corresponding to internal material numbers
         return true;
      }
      //--------------------------------------------------------------------
      test="minimum-height";
      if(word==test){
          double min=atof(value.c_str());
          vin::check_for_valid_value(min, word, line, prefix, unit, "none", 0.0, 1.0,"material"," 0.0 - 1.0");
          create::internal::select_material_by_z_height = true; // default
          create::internal::mp[super_index].min=min;
          return true;
      }
      //--------------------------------------------------------------------
      test="maximum-height";
      if(word==test){
          double max=atof(value.c_str());
          vin::check_for_valid_value(max, word, line, prefix, unit, "none", 0.0, 1.0,"material"," 0.0 - 1.0");
          create::internal::select_material_by_z_height = true; // default
          create::internal::mp[super_index].max=max;
          return true;
      }
      //--------------------------------------------------------------------
      // keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of namespace create
