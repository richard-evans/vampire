//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <sstream>

// Vampire headers
#include "errors.hpp"
#include "sim.hpp"
#include "vio.hpp"

// Internal sim header
#include "internal.hpp"

namespace sim{

   //-----------------------------------------------------------------------------
   // Function to process input file parameters for sim module
   //-----------------------------------------------------------------------------
   bool match_input_parameter(string const key, string const word, string const value, string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="sim";
      if(key!=prefix) return false;

      // set maximum allowable value for time steps (10^12)
      const uint64_t max_time = 1000000000000;
      const std::string max_time_str = "1,000,000,000,000";

      //----------------------------------
      // Now test for all valid options
      //----------------------------------
      std::string test="slonczewski-spin-polarization-unit-vector";
      if(word==test){
         std::vector<double> u(3);
         u=vin::doubles_from_string(value);
         // Test for valid range
         vin::check_for_valid_unit_vector(u, word, line, prefix, "input");
         // save sanitized unit vector
         sim::internal::slonczewski_spin_polarization_unit_vector = u;
         return true;
      }
      //-------------------------------------------------------------------
      test="SOT-spin-polarization-unit-vector";
      if(word==test){
         std::vector<double> u(3);
         u=vin::doubles_from_string(value);
         // Test for valid range
         vin::check_for_valid_unit_vector(u, word, line, prefix, "input");
         // save sanitized unit vector
         sim::internal::SOT_spin_polarization_unit_vector = u;
         return true;
      }
      //-------------------------------------------------------------------
      test="preconditioning-steps";
      if(word==test){
         int n = atoi(value.c_str());
         // Test for valid range
         vin::check_for_valid_int(n, word, line, prefix, 0, 1000000,"input","0 - 1,000,000");
         sim::num_monte_carlo_preconditioning_steps = n;
         return true;
      }
      //-------------------------------------------------------------------
      test="time-step";
      if(word==test){
         double dt = atof(value.c_str());
         vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-6,"input","0.01 attosecond - 1 picosecond");
         mp::dt_SI = dt;
         return true;
      }
      //--------------------------------------------------------------------
      test="total-time-steps";
      if(word==test){
         uint64_t tt = vin::str_to_uint64(value); // convert string to uint64_t
         vin::check_for_valid_int(tt, word, line, prefix, 0, max_time,"input","0 - "+max_time_str);
         sim::total_time = tt;
         return true;
      }
      //--------------------------------------------------------------------
      test="loop-time-steps";
      if(word==test){
         uint64_t tt = vin::str_to_uint64(value); // convert string to uint64_t
         vin::check_for_valid_int(tt, word, line, prefix, 0, max_time,"input","0 - "+max_time_str);
         sim::loop_time = tt;
         return true;
      }
      //--------------------------------------------------------------------
      test="equilibration-time-steps";
      if(word==test){
         uint64_t tt = vin::str_to_uint64(value); // convert string to uint64_t
         vin::check_for_valid_int(tt, word, line, prefix, 0, max_time,"input","0 - "+max_time_str);
         sim::equilibration_time = tt;
         return true;
      }
      //--------------------------------------------------------------------
      test="time-steps-increment";
      if(word==test){
         uint64_t tt = vin::str_to_uint64(value); // convert string to uint64_t
         vin::check_for_valid_int(tt, word, line, prefix, 1, max_time,"input","1 - "+max_time_str);
         sim::partial_time = tt;
         return true;
      }
      test="domain-wall-axis";
      if(word==test){
         //vin::check_for_valid_int(tt, word, line, prefix, 0, max_time,"input","0 - "+max_time_str);
         if (value == "x") {
         sim::domain_wall_axis = 0;
         }
         else if (value == "y") {
         sim::domain_wall_axis = 1;
         }
         else if (value == "z") {
         sim::domain_wall_axis = 2;
         }
         else {
            std::cout << "domain wall axis must equal x or y or z" <<std::endl;
            return false;
         }
         return true;
      }
      test="domain-wall-discretisation";
      if(word==test){
         double tt = atof(value.c_str()); // convert string to uint64_t
         vin::check_for_valid_value(tt, word, line, prefix, unit, "length", 10, 1000,"input","10 - 1 A");
         sim::domain_wall_discretisation = tt;
         return true;
      }
      test="domain-wall-anti-pbc-x";
      if(word==test){
         sim::anti_PBC[0] = true;
         cs::pbc[0]=true;
         return true;
      }
      test="domain-wall-anti-pbc-y";
      if(word==test){
         sim::anti_PBC[1] = true;
         cs::pbc[1]=true;
         return true;
      }
      test="domain-wall-anti-pbc-z";
      if(word==test){
         sim::anti_PBC[2] = true;
         cs::pbc[2]=true;
         return true;
      }

      test="domain-wall-position";
      if(word==test){
         double tt = atof(value.c_str()); // convert string to uint64_t
         vin::check_for_valid_value(tt, word, line, prefix, unit, "none", 0, 1,"input","0 - 1");
         sim::domain_wall_position = tt;
         return true;
      }
      test="domain-wall-width";
      if(word==test){
         double tt = atof(value.c_str()); // convert string to uint64_t
         vin::check_for_valid_value(tt, word, line, prefix, unit, "length", 0, 1000,"input","0 - 1 A");
         sim::domain_wall_width = tt;
         return true;
      }


      //--------------------------------------------------------------------
      // input parameter not found here
      return false;
   }

   //----------------------------------------------------------------------------------
   // material parameter match function
   //----------------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index){

      // add prefix string
      std::string prefix="material:";

      // Check for material id > current array size and if so dynamically expand mp array
      if((unsigned int) super_index + 1 > sim::internal::mp.size() && super_index + 1 < 101) sim::internal::mp.resize(super_index + 1);

      //------------------------------------------------------------
      std::string test="slonczewski-adiabatic-spin-torque";
      /*
         aj parameter for material in slonczewski torque calculation
         */
      if(word==test){
         double aj=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(aj, word, line, prefix, unit, "field", 0.0, 1.0e2,"input","0 - 100T");
         sim::internal::mp[super_index].slonczewski_aj.set(aj);
         return true;
      }
      test="domain-wall-second-magnetisation-vector";
      if(word==test){
       std::vector<double> u(3);
       u=vin::doubles_from_string(value);
        vin::check_for_valid_unit_vector(u, word, line, prefix, "input");
        std::cout << sim::domain_wall_second_vector_x.size() << "\t" << super_index << "\t" << u[0] << '\t' << u[1] << '\t' << u[2] <<std::endl;
        sim::domain_wall_second_vector_x[super_index] = u[0];
        sim::domain_wall_second_vector_y[super_index] = u[1];
        sim::domain_wall_second_vector_z[super_index] = u[2];
         return true;
      }
      //------------------------------------------------------------
      test="slonczewski-non-adiabatic-spin-torque";
      /*
         bj parameter for material in slonczewski torque calculation
         */
      if(word==test){
         double bj=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(bj, word, line, prefix, unit, "field", 0.0, 1.0e2,"input","0 - 100T");
         sim::internal::mp[super_index].slonczewski_bj.set(bj);
         return true;
      }
      //------------------------------------------------------------
      test="SOT-damping-like-torque";
      /*
         damping-like parameter for material in spin orbit torque calculation
         */
      if(word==test){
         double aj=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(aj, word, line, prefix, unit, "field", -1.0e2, 1.0e2,"input","-100 - 100T");
         sim::internal::mp[super_index].SOT_DL.set(aj);
         return true;
      }
      //------------------------------------------------------------
      test="SOT-field-like-torque";
      /*
         field-like parameter for material in spin orbit torque calculation
         */
      if(word==test){
         double bj=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(bj, word, line, prefix, unit, "field", -1.0e2, 1.0e2,"input","-100 - 100T");
         sim::internal::mp[super_index].SOT_FL.set(bj);
         return true;
      }



      //--------------------------------------------------------------------
      // keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of namespace sim
