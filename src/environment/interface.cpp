//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins 2016. All rights reserved.
//
//   Email: sj681@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <string>

// Vampire headers
#include "environment.hpp"
#include "errors.hpp"
#include "vio.hpp"

// environment module headers
#include "internal.hpp"

namespace env = environment::internal;

namespace environment{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for environment module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="environment";
      if(key!=prefix) return false;

      std::string test="enabled";
      if(word==test){
        enabled = true;
        return true;
      }

      test="system-size-x";
      if(word==test){
        double x=atof(value.c_str());
        //vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-6,"input","0.01 attosecond - 1 picosecond");
        env::dim[0] =x;
        return true;
      }
      test="system-size-y";
      if(word==test){
        double y=atof(value.c_str());
        //vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-6,"input","0.01 attosecond - 1 picosecond");
        env::dim[1] =y;
        return true;
      }
      test="system-size-z";
      if(word==test){
        double z=atof(value.c_str());
        //vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-6,"input","0.01 attosecond - 1 picosecond");
        env::dim[2] =z;
        return true;
      }
      test="system-size";
      if(word==test){
        double d=atof(value.c_str());
        //vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-6,"input","0.01 attosecond - 1 picosecond");
        env::dim[0] =d;
        env::dim[1] =d;
        env::dim[2] =d;
        return true;
      }

      test="cell-size";
      if(word==test){
        double d=atof(value.c_str());
        //vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-6,"input","0.01 attosecond - 1 picosecond");
        env::cell_size[0] =d;
        env::cell_size[1] =d;
        env::cell_size[2] =d;
        return true;
      }
      test="cell-size-x";
      if(word==test){
        double d=atof(value.c_str());
        //vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-6,"input","0.01 attosecond - 1 picosecond");
        env::cell_size[0] =d;
        return true;
      }
      test="cell-size-y";
      if(word==test){
        double d=atof(value.c_str());
        //vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-6,"input","0.01 attosecond - 1 picosecond");
        env::cell_size[1] =d;
        return true;
      }
      test="cell-size-z";
      if(word==test){
        double d=atof(value.c_str());
        //vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-6,"input","0.01 attosecond - 1 picosecond");
        env::cell_size[2] =d;
        return true;
      }

      test="exchange-constant";
      if(word==test){
        double exc=atof(value.c_str());
        //vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-6,"input","0.01 attosecond - 1 picosecond");
        env::A =-exc;
        return true;
      }

      test="uniaxial-anisotropy-constant";
      if(word==test){
        double k=atof(value.c_str());
        //vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-6,"input","0.01 attosecond - 1 picosecond");
        env::ku =k;
        return true;
      }
      test="Ms";
      if(word==test){
        double A=atof(value.c_str());
        //vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-6,"input","0.01 attosecond - 1 picosecond");
        env::Ms =A;
        return true;
      }


      test="Tc";
      if(word==test){
        double T=atof(value.c_str());
        //vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-6,"input","0.01 attosecond - 1 picosecond");
        env::Tc =T;
        return true;
      }

      test="damping-constant";
      if(word==test){
        double a=atof(value.c_str());
        //vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-6,"input","0.01 attosecond - 1 picosecond");
        env::alpha =a;
        return true;
      }
      test="gyromagnetic-ratio";
      if(word==test){
        double g=atof(value.c_str());
        //vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-6,"input","0.01 attosecond - 1 picosecond");
        env::gamma =g;
        return true;
      }
      test="initial-spin-direction";
      if(word==test){
          // first test for random spins
          test="random";
          if(value==test){
              env::random_spins=true;
          }
          else{
              // temporary storage container
              std::vector<double> u(3);

              // read values from string
              u=vin::DoublesFromString(value);


              // Copy sanitised unit vector to material
              env::initial_spin[0]=u.at(0);
              env::initial_spin[1]=u.at(1);
              env::initial_spin[2]=u.at(2);

              // ensure random spins is unset
              env::random_spins=false;
          }
          // return
          return true;
      }

      test="atomistic-position";
      if(word==test){
        // temporary storage container
        std::vector<double> u(3);

        // read values from string
        u=vin::DoublesFromString(value);
        // Copy sanitised unit vector to material
        env::shift[0]=u.at(0);
        env::shift[1]=u.at(1);
        env::shift[2]=u.at(2);
        return true;
      }
      test="demag-field-update-rate";
      if(word==test){
         int dpur=atoi(value.c_str());
         vin::check_for_valid_int(dpur, word, line, prefix, 0, 1000000,"input","0 - 1,000,000");
         environment::demag_update_rate=dpur;
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

} // end of environment namespace
