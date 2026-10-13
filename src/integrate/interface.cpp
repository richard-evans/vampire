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

// Vampire headers
#include "integrate.hpp"
#include "errors.hpp"
#include "random.hpp"
#include "vio.hpp"

// integrate module headers
#include "internal.hpp"

namespace integrate{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for integrate module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      //------------------------------------------------------------------------
      // Integrator keywords under the sim: prefix (kept for compatibility
      // with existing input files)
      //------------------------------------------------------------------------
      if(key == "sim"){

         std::string prefix="sim";
         std::string test;

         test = "integrator";
         if( word == test ){
            //--------------------------------------------------------------------
            test="llg-heun";
            if( value == test ){
               integrate::integrator = integrator_t::llg_heun;
               return true;
            }
            //--------------------------------------------------------------------
            test="monte-carlo";
            if( value == test ){
               integrate::integrator = integrator_t::monte_carlo;
               return true;
            }
            //--------------------------------------------------------------------
            test="llg-midpoint";
            if( value == test ){
               integrate::integrator = integrator_t::llg_midpoint;
               return true;
            }
            //--------------------------------------------------------------------
            test="constrained-monte-carlo";
            if( value == test ){
               integrate::integrator = integrator_t::cmc;
               return true;
            }
            //--------------------------------------------------------------------
            test="hybrid-constrained-monte-carlo";
            if( value == test ){
               integrate::integrator = integrator_t::hybrid_cmc;
               return true;
            }
            //--------------------------------------------------------------------
            test="llg-quantum";
            if( value == test ){
               integrate::integrator = integrator_t::llg_quantum;
               return true;
            }
            //--------------------------------------------------------------------
            test="llg-rkmk2";
            if( value == test ){
               integrate::integrator = integrator_t::llg_rkmk2;
               return true;
            }
            //--------------------------------------------------------------------
            test="llg-rkmk4";
            if( value == test ){
               integrate::integrator = integrator_t::llg_rkmk4;
               return true;
            }
            //--------------------------------------------------------------------
            test="llg-dm";
            if( value == test ){
               integrate::integrator = integrator_t::llg_dm;
               return true;
            }
            //--------------------------------------------------------------------
            test="spin-lattice";
            if( value == test ){
               integrate::integrator = integrator_t::suzuki_trotter;
               return true;
            }
            //
            //--------------------------------------------------------------------
            test="lsf";
            if( value == test ){
               integrate::integrator = integrator_t::lsf;
               return true;
            }
            //--------------------------------------------------------------------
            test="lsf-mc";
            if( value == test ){
               integrate::integrator = integrator_t::lsf_mc;
               return true;
            }
            //--------------------------------------------------------------------
            test="lsf-rk4";
            if( value == test ){
               integrate::integrator = integrator_t::lsf_rk4;
               return true;
            }
            //--------------------------------------------------------------------
            else{
               terminaltextcolor(RED);
                  std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
                  std::cerr << "\t\"llg-heun\"" << std::endl;
                  std::cerr << "\t\"llg-midpoint\"" << std::endl;
                  std::cerr << "\t\"llg-quantum\"" << std::endl;
                  std::cerr << "\t\"llg-rkmk2\"" << std::endl;
                  std::cerr << "\t\"llg-rkmk4\"" << std::endl;
                  std::cerr << "\t\"llg-dm\"" << std::endl;
                  std::cerr << "\t\"monte-carlo\"" << std::endl;
                  std::cerr << "\t\"constrained-monte-carlo\"" << std::endl;
                  std::cerr << "\t\"lsf\"" << std::endl;
                  std::cerr << "\t\"lsf-mc\"" << std::endl;
                  std::cerr << "\t\"lsf-rk4\"" << std::endl;
               terminaltextcolor(WHITE);
               err::vexit();
             }
         }
         //--------------------------------------------------------------------
         test="integrator-random-seed";
         if(word==test){
            int is=atoi(value.c_str());
            vin::check_for_valid_int(is, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
            mtrandom::integration_seed=is;
            return true;
         }
         //--------------------------------------------------------------------
         // sim: keyword not found here
         //--------------------------------------------------------------------
         return false;

      }

      // Check for valid key, if no match return false
      std::string prefix="integrate";
      if(key!=prefix) return false;

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

      // Check for material id > current array size and if so dynamically expand mp array
      if((unsigned int) super_index + 1 > internal::mp.size() && super_index + 1 < 101) internal::mp.resize(super_index + 1);

      std::string test;

      //------------------------------------------------------------
      test="lsf-second-order-coefficient"; // Landau coefficient A
      if( word==test ){
         double a_l=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(a_l, word, line, prefix, unit, "energy", -1.0e-16, 1.0e-16,"material","-1.0e-16 - 1.0e-16");
         internal::mp[super_index].lsf_second_order_coefficient.set(a_l);
         return true;
      }
      //------------------------------------------------------------
      test="lsf-fourth-order-coefficient"; // Landau coefficient B
      if( word==test ){
         double b_l=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(b_l, word, line, prefix, unit, "energy", -1.0e-16, 1.0e-16,"material","-1.0e-16 - 1.0e-16");
         internal::mp[super_index].lsf_fourth_order_coefficient.set(b_l);
         return true;
      }
      //------------------------------------------------------------
      test="lsf-sixth-order-coefficient"; // Landau coefficient C
      if( word==test ){
         double c_l=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(c_l, word, line, prefix, unit, "energy", -1.0e-16, 1.0e-16,"material","-1.0e-16 - 1.0e-16");
         internal::mp[super_index].lsf_sixth_order_coefficient.set(c_l);
         return true;
      }

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of integrate namespace

