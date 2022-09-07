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
#include "voronoi.hpp"
#include "random.hpp"
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
         u=vin::doubles_from_string(value);
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

      // cs::system_creation_flags needs refactoring for readability and bug resistance
      test="full";
      if(word==test){
         cs::system_creation_flags[1]=0;
         return true;
      }
      //-------------------------------------------------------------------
      test="cube";
      if(word==test){
         cs::system_creation_flags[1]=1;
         return true;
      }
      //-------------------------------------------------------------------
      test="cylinder";
      if(word==test){
         cs::system_creation_flags[1]=2;
         return true;
      }
      //-------------------------------------------------------------------
      test="ellipsoid";
      if(word==test){
         cs::system_creation_flags[1]=3;
         return true;
      }
      //-------------------------------------------------------------------
      test="sphere";
      if(word==test){
         cs::system_creation_flags[1]=4;
         return true;
      }
      //-------------------------------------------------------------------
      test="truncated-octahedron";
      if(word==test){
         cs::system_creation_flags[1]=5;
         return true;
      }
      //-------------------------------------------------------------------
      test="tear-drop";
      if(word==test){
         cs::system_creation_flags[1]=6;
         return true;
      }
      //-------------------------------------------------------------------
      // system_creation_flags[2] - Set system type
      //-------------------------------------------------------------------
      test="particle";
      if(word==test){
         cs::system_creation_flags[2]=0;
         return true;
      }
      else
      test="particle-array";
      if(word==test){
         cs::system_creation_flags[2]=1;
         return true;
      }
      else
      test="hexagonal-particle-array";
      if(word==test){
         cs::system_creation_flags[2]=2;
         return true;
      }
      else
      test="voronoi-film";
      if(word==test || word == "grains" || word == "granular-film"){
         cs::system_creation_flags[2]=3;
         return true;
      }
      test="voronoi-grain-substructure";
      if(word==test){
         create::internal::generate_voronoi_substructure = true;
         return true;
      }
      ///-------------------------------------------------------------------
      /// system_creation_flags[1] - Set system particle shape
      ///-------------------------------------------------------------------

      //--------------------------------------------------------------------
      test="voronoi-size-variance";
      if(word==test || word == "grain-size-variance"){
         double vsd=atof(value.c_str());
         vin::check_for_valid_value(vsd, word, line, prefix, unit, "none", 0.0, 1.0,"input","0.0 - 1.0");
         create_voronoi::voronoi_sd=vsd;
         return true;
      }
      //--------------------------------------------------------------------
      test="voronoi-row-offset";
      if(word==test){
         create_voronoi::parity=1;
         return true;
      }
      //--------------------------------------------------------------------
      test="voronoi-random-seed";
      if(word==test){
         int vs=atoi(value.c_str());
         vin::check_for_valid_int(vs, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
         mtrandom::voronoi_seed=vs;
         return true;
      }
      test="voronoi-rounded-grains";
      if(word==test || word == "rounded-grains"){
         create_voronoi::rounded=true;
         return true;
      }
      test="voronoi-include-boundary-grains";
      if(word==test || word == "include-boundary-grains"){
         create_voronoi::include_boundary_grains_real=true;
         return true;
      }
      //-------------------------------------------------------------------
      test="voronoi-rounded-grains-area";
      if(word==test || word == "rounded-grains-area"){
         double vsd=atof(value.c_str());
         vin::check_for_valid_value(vsd, word, line, prefix, unit, "none", 0.0, 1.0,"input","0.0 - 1.0");
         create_voronoi::area_cutoff=vsd;
         return true;
      }
      //-------------------------------------------------------------------
      test="particle-centre-offset"; //parity
      if(word==test){
         cs::particle_creation_parity=1;
         return true;
      }
      //--------------------------------------------------------------------
      test="single-spin";
      if(word==test){
         cs::single_spin=true;
         return true;
      }
      //--------------------------------------------------------------------
      test="periodic-boundaries";
      if( word == test ){
         // set PBC for all directions
         if(value == ""){
            cs::pbc[0] = true;
            cs::pbc[1] = true;
            cs::pbc[2] = true;
         }
         // search for x,y,z in value
         std::string strs[3] = { "x", "y", "z" };
         for(int i=0; i<3; i++){
            if (value.find(strs[i]) != std::string::npos) cs::pbc[i] = true;
         }
         // check for unknown junk in value
         if( value.find_first_not_of("xyz,") != std::string::npos){
            terminaltextcolor(RED);
            std::cerr << "Error - value for \'create:" << word << "\' must be a combination of x,y,z only" << std::endl;
            zlog << zTs() << "Error - value for \'create:" << word << "\' must be a combination of x,y,z only" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
            return false;
         }
         return true;
      }
      //--------------------------------------------------------------------
      test="periodic-boundaries-x";
      if(word==test){
         cs::pbc[0]=true;
         return true;
      }
      //--------------------------------------------------------------------
      test="periodic-boundaries-y";
      if(word==test){
         cs::pbc[1]=true;
         return true;
      }
      //--------------------------------------------------------------------
      test="periodic-boundaries-z";
      if(word==test){
         cs::pbc[2]=true;
         return true;
      }
      //--------------------------------------------------------------------
      test="select-material-by-geometry";
      if(word==test){
         create::internal::select_material_by_geometry = true; // default
         // also check for value
         std::string vfalse="false";
         if(value == vfalse){
            create::internal::select_material_by_geometry = false;
         }
         return true;
      }
      //--------------------------------------------------------------------
      test="fill-core-shell-particles";
      if(word==test){
         cs::fill_core_shell=true; // default
         // also check for value
         std::string VFalse="false";
         if(value==VFalse){
            cs::fill_core_shell=false;
         }
         return true;
      }
      //--------------------------------------------------------------------
      test="interfacial-roughness";
      if(word==test){
         cs::interfacial_roughness=true; // default
         // also check for value
         std::string VFalse="false";
         if(value==VFalse){
            cs::interfacial_roughness=false;
         }
         return true;
      }
      //--------------------------------------------------------------------
      test="material-interfacial-roughness";
      if(word==test){
         cs::interfacial_roughness_local_height_field=true; // default
         // also check for value
         std::string VFalse="false";
         if(value==VFalse){
            cs::interfacial_roughness_local_height_field=false;
         }
         return true;
      }
      //--------------------------------------------------------------------
      test="interfacial-roughness-random-seed";
      if(word==test){
         unsigned int vs=atoi(value.c_str());
         cs::interfacial_roughness_random_seed=vs;
         return true;
      }
      //--------------------------------------------------------------------
      test="interfacial-roughness-number-of-seed-points";
      if(word==test){
         int sc=atoi(value.c_str());
         vin::check_for_valid_int(sc, word, line, prefix, 0, 100000,"input","0 - 100,000");
         cs::interfacial_roughness_seed_count=sc;
         return true;
      }
      //--------------------------------------------------------------------
      test="interfacial-roughness-type";
      if(word==test){
         std::string loctest="peaks";
         if(value==loctest){
            cs::interfacial_roughness_type=1;
            return true;
         }
         else
         loctest="troughs";
         if(value==loctest){
            cs::interfacial_roughness_type=-1;
            return true;
         }
         else{
            cs::interfacial_roughness_type=0;
            return true;
         }
      }
      //--------------------------------------------------------------------
      test="interfacial-roughness-seed-radius";
      if(word==test){
         double irsr=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(irsr, word, line, prefix, unit, "length", 0.0, 10000.0,"input","0.0 - 1 micrometre");
         cs::interfacial_roughness_mean_seed_radius=irsr;
         return true;
      }
      //--------------------------------------------------------------------
      test="interfacial-roughness-seed-radius-variance";
      if(word==test){
         double irsrv=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(irsrv, word, line, prefix, unit, "none", 0.0, 1.0,"input","0.0 - 1.0");
         cs::interfacial_roughness_seed_radius_variance=irsrv;
         return true;
      }
      //--------------------------------------------------------------------
      test="interfacial-roughness-mean-height";
      if(word==test){
         double irmh=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(irmh, word, line, prefix, unit, "length", 0.1, 100.0,"input","0.1 Angstroms - 10 nanometres");
         cs::interfacial_roughness_mean_seed_height=irmh;
         return true;
      }
      //--------------------------------------------------------------------
      test="interfacial-roughness-maximum-height";
      if(word==test){
         double shm=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(shm, word, line, prefix, unit, "length", 0.1, 100.0,"input","0.1 Angstroms - 10 nanometres");
         cs::interfacial_roughness_seed_height_max=shm;
         return true;
      }
      //--------------------------------------------------------------------
      test="interfacial-roughness-height-field-resolution";
      if(word==test){
         double irhfr=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(irhfr, word, line, prefix, unit, "length", 0.1, 100.0,"input","0.1 Angstroms - 10 nanometres");
         cs::interfacial_roughness_height_field_resolution=irhfr;
         return true;
      }
      //--------------------------------------------------------------------
      test="multilayers";
      if(word==test){
         int nmul=atoi(value.c_str());
         // Test for valid range
         vin::check_for_valid_int(nmul, word, line, prefix, 1, 100,"input","1 - 100, specifying the number of multilayers to be generated");
         cs::multilayers = true;
         cs::num_multilayers = nmul;
         return true;
      }
      //--------------------------------------------------------------------
      test="grains-poission"; //parity
      if(word==test){
         create::internal::grain_poission=true;
         return true;
      }
      //--------------------------------------------------------------------
      test="height-categorization";
      if(word==test){
         // Test for different options
         test="default";
         if(value==test){
            // do nothing
            return true;
         }
         test="multilayers";
         if(value==test){
            cs::multilayer_height_category = true;
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error - value for \'create:" << word << "\' must be one of:" << std::endl;
            std::cerr << "\t\"default\"" << std::endl;
            std::cerr << "\t\"multilayers\"" << std::endl;
            zlog << zTs() << "Error - value for \'create:" << word << "\' must be one of:" << std::endl;
            zlog << zTs() << "\t\"default\"" << std::endl;
            zlog << zTs() << "\t\"multilayers\"" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }
      //--------------------------------------------------------------------
      test="voronoi-grain-substructure-crystallization-radius";
      if(word==test){
         double rsize=atof(value.c_str());
         vin::check_for_valid_value(rsize, word, line, prefix, unit, "none", 0.01, 2.0,"input","0.01 - 2");
         create::internal::voronoi_grain_substructure_crystallization_radius=rsize;
         return true;
      }
      //--------------------------------------------------------------------
      test="voronoi-grain-substructure-overlap-factor";
      if(word==test){
         double ol=atof(value.c_str());
         vin::check_for_valid_value(ol, word, line, prefix, unit, "none", 0.1, 3.0,"input","0.1 - 3");
         create::internal::voronoi_grain_substructure_overlap_factor = ol;
         return true;
      }
      //--------------------------------------------------------------------
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
      test="ellipse";
      if(word==test){
         cs::system_creation_flags[1] = 10;
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
         return true;
      }
      //--------------------------------------------------------------------
      test="alloy-random-seed";
      if(word==test){
         int ars=atoi(value.c_str());
         vin::check_for_valid_int(ars, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
         create::internal::alloy_seed = ars;
         return true;
      }
      //--------------------------------------------------------------------
      test="grain-random-seed";
      if(word==test){
         int grs=atoi(value.c_str());
         vin::check_for_valid_int(grs, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
         create::internal::grain_seed = grs;
         return true;
      }
      //--------------------------------------------------------------------
      test="dilution-random-seed";
      if(word==test){
         int drs=atoi(value.c_str());
         vin::check_for_valid_int(drs, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
         create::internal::dilute_seed = drs;
         return true;
      }
      //--------------------------------------------------------------------
      test="intermixing-random-seed";
      if(word==test){
         int mrs=atoi(value.c_str());
         vin::check_for_valid_int(mrs, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
         create::internal::mixing_seed = mrs;
         return true;
      }
      //--------------------------------------------------------------------
      test="spin-initialisation-random-seed";
      if(word==test){
         int sirs=atoi(value.c_str());
         vin::check_for_valid_int(sirs, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
         create::internal::spin_init_seed = sirs;
         return true;
      }
      /*std::string test="slonczewski-spin-polarization-unit-vector";
      if(word==test){
         std::vector<double> u(3);
         u=vin::doubles_from_string(value);
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
      //--------------------------------------------------------------------
      test="geometry-file";
      if( word == test ){

         //------------------------------------------------------------------
         // Open geometry file
         //------------------------------------------------------------------
         std::stringstream gfile;
         gfile.str( vin::get_string(value.c_str(), "material", line) );

         //------------------------------------------------------------------
         // read in number of points in geometry file
         //------------------------------------------------------------------
         int num_points = 0; // number of points in geometry file
         gfile >> num_points;


         //------------------------------------------------------------------
         // check for a sensible number of points
         //------------------------------------------------------------------
         if( (num_points < 3) || (num_points > 1000) ){
            std::cerr << "Error in geometry input file " << value.c_str() << " - first number must be non zero integer in the range 3-100"<< std::endl;
            err::vexit();
         }

         //------------------------------------------------------------------
         // reserve memory for points
         //------------------------------------------------------------------
         create::internal::mp[super_index].geometry_points.reserve(num_points);

         //------------------------------------------------------------------
         // Process points in file
         //------------------------------------------------------------------
         for( int c = 0; c < num_points; c++){

            double pts[2] = { 0.0, 0.0 }; // temporary variables for points
            // read one coordinate at a time
            for( int xy = 0 ; xy < 2 ; xy++ ){
               double var;
               gfile >> var;
               // check for no end of file
               if(gfile.eof()){
                  std::cerr << "Error in geometry input file " << value.c_str() << " end of file reached before reading all coordinates" << std::endl;
                  err::vexit();
               }
               // check for sensible coordinate
               if((var<0.0) || (var > 1.0)){
                  std::cerr << "Error in geometry input file " << value.c_str() << " value is outside of valid range (0.0-1.0)" << std::endl;
                  err::vexit();
               }
               // now save sanitised input
               pts[xy] = var;
            }
            // save points in materials class using constructor
            create::internal::mp[super_index].geometry_points.emplace_back(pts[0],pts[1]);

         }

         // enable geometry
         create::internal::mp[super_index].geometry = true;

         return true;

      }
      //--------------------------------------------------------------------
      test="voronoi-grain-substructure-nucleation-height";
      if(word==test){
         double nh=atof(value.c_str());
         vin::check_for_valid_value(nh, word, line, prefix, unit, "none", 0.0, 1.0,"material"," 0.0 - 1.0");
         create::internal::mp[super_index].voronoi_grain_substructure_nucleation_height = nh;
         return true;
      }
      //--------------------------------------------------------------------
      /*
         integer to associate the material to a particular material within the unit cell.
         Default is 0 but can be overidden with this parameter.
      */
      test="unit-cell-category";
      if(word==test){
         int uccat=atoi(value.c_str());
         vin::check_for_valid_int(uccat, word, line, prefix, 1, mp::max_materials+1,"material"," 1 - 100");
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
