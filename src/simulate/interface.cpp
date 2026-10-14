//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2025. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <fstream>
#include <sstream>

// Vampire headers
#include "errors.hpp"
#include "integrate.hpp"
#include "montecarlo.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "units.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// Internal sim header
#include "internal.hpp"
#include "sld.hpp"

namespace sim{

   //-----------------------------------------------------------------------------
   // Function to process input file parameters for sim module
   //-----------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="sim";
      if(key!=prefix) return false;

      // set maximum allowable value for time steps (10^12)
      const uint64_t max_time = 1000000000000;
      const std::string max_time_str = "1,000,000,000,000";

      //----------------------------------
      // Now test for all valid options
      //----------------------------------
      std::string test="spin-transfer-torque-polarization-unit-vector";
      if(word==test){
         std::vector<double> u(3);
         u=vin::doubles_from_string(value);
         // Test for valid range
         vin::check_for_valid_unit_vector(u, word, line, prefix, "input");
         // save sanitized unit vector
         sim::internal::stt_polarization_unit_vector.x = u[0];
         sim::internal::stt_polarization_unit_vector.y = u[1];
         sim::internal::stt_polarization_unit_vector.z = u[2];
         return true;
      }
      //-------------------------------------------------------------------
      test="spin-orbit-torque-polarization-unit-vector";
      if(word==test){
         std::vector<double> u(3);
         u=vin::doubles_from_string(value);
         // Test for valid range
         vin::check_for_valid_unit_vector(u, word, line, prefix, "input");
         // save sanitized unit vector
         sim::internal::sot_polarization_unit_vector.x = u[0];
         sim::internal::sot_polarization_unit_vector.y = u[1];
         sim::internal::sot_polarization_unit_vector.z = u[2];
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
      //--------------------------------------------------------------------
      test="domain-wall-axis";
      if(word==test){
         //vin::check_for_valid_int(tt, word, line, prefix, 0, max_time,"input","0 - "+max_time_str);
         if (value == "x") {
         program::domain_wall_axis = 0;
         }
         else if (value == "y") {
         program::domain_wall_axis = 1;
         }
         else if (value == "z") {
         program::domain_wall_axis = 2;
         }
         else {
            std::cout << "domain wall axis must equal x or y or z" <<std::endl;
            return false;
         }
         return true;
      }
      //--------------------------------------------------------------------
      test="domain-wall-discretisation";
      if(word==test){
         double tt = atof(value.c_str()); // convert string to uint64_t
         vin::check_for_valid_value(tt, word, line, prefix, unit, "length", 10, 1000,"input","10 - 1 A");
         program::domain_wall_discretisation = tt;
         return true;
      }
      //--------------------------------------------------------------------
      test="domain-wall-anti-pbc-x";
      if(word==test){
         program::anti_PBC[0] = true;
         cs::pbc[0]=true;
         return true;
      }
      //--------------------------------------------------------------------
      test="domain-wall-anti-pbc-y";
      if(word==test){
         program::anti_PBC[1] = true;
         cs::pbc[1]=true;
         return true;
      }
      //--------------------------------------------------------------------
      test="domain-wall-anti-pbc-z";
      if(word==test){
         program::anti_PBC[2] = true;
         cs::pbc[2]=true;
         return true;
      }
      //--------------------------------------------------------------------
      test="domain-wall-position";
      if(word==test){
         double tt = atof(value.c_str()); // convert string to uint64_t
         vin::check_for_valid_value(tt, word, line, prefix, unit, "none", 0, 1,"input","0 - 1");
         program::domain_wall_position = tt;
         return true;
      }
      //--------------------------------------------------------------------
      test="domain-wall-width";
      if(word==test){
         double tt = atof(value.c_str()); // convert string to uint64_t
         vin::check_for_valid_value(tt, word, line, prefix, unit, "length", 0, 1000,"input","0 - 1 A");
         program::domain_wall_width = tt;
         return true;
      }
      //--------------------------------------------------------------------
      test="enable-fmr-field";
      if(word==test){
         return true;
      }
      //--------------------------------------------------------------------
      test="simulation-cycles";
      if(word==test){
         int r=atoi(value.c_str());
         vin::check_for_valid_int(r, word, line, prefix, 0, 2000000000,"input","1 - 2,000,000,000");
         sim::runs=r;
         return true;
      }
      //--------------------------------------------------------------------
      test="maximum-temperature";
      if(word==test){
         double T=atof(value.c_str());
         vin::check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 1.0e6,"input","0.0 - 1,000,000 K");
         sim::Tmax=T;
         return true;
      }
      //--------------------------------------------------------------------
      test="minimum-temperature";
      if(word==test){
         double T=atof(value.c_str());
         vin::check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 1.0e6,"input","0.0 - 1,000,000 K");
         sim::Tmin=T;
         return true;
      }
      //--------------------------------------------------------------------
      test="equilibration-temperature";
      if(word==test){
         double T=atof(value.c_str());
         vin::check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 1.0e6,"input","0.0 - 1,000,000 K");
         sim::Teq=T;
         return true;
      }
      //--------------------------------------------------------------------
      test="temperature";
      if(word==test){
         double T=atof(value.c_str());
         vin::check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 1.0e6,"input","0.0 - 1,000,000 K");
         sim::temperature=T;
         return true;
      }
      //--------------------------------------------------------------------
      test="temperature-increment";
      if(word==test){
         double T=atof(value.c_str());
         vin::check_for_valid_value(T, word, line, prefix, unit, "none", 1.0e-10, 1.0e6,"input","1.0e-10 - 1,000,000 K");
         sim::delta_temperature=T;
         return true;
      }
      //--------------------------------------------------------------------
      test="cooling-time";
      if(word==test){
         double T=atof(value.c_str());
         vin::check_for_valid_value(T, word, line, prefix, unit, "time", 1.0e-18, 1.0,"input","1 attosecond - 1 s");
         program::cooling_time=T;
         return true;
      }
      //--------------------------------------------------------------------
      test="cooling-function";
      if(word==test){
         test="exponential";
         if(value==test){
            program::cooling_function_flag=0;
            return true;
         }
         test="gaussian";
         if(value==test){
            program::cooling_function_flag=1;
            return true;
         }
         test="double-gaussian";
         if(value==test){
            program::cooling_function_flag=2;
            return true;
         }
         test="linear";
         if(value==test){
            program::cooling_function_flag=3;
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
            std::cerr << "\t\"exponential\"" << std::endl;
            std::cerr << "\t\"gaussian\"" << std::endl;
            std::cerr << "\t\"double-gaussian\"" << std::endl;
            std::cerr << "\t\"linear\"" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }
      //--------------------------------------------------------------------
      test="laser-pulse-temporal-profile";
      if(word==test){
         test="square";
         if(value==test){
            program::pump_function=program::square;
            return true;
         }
         test="two-temperature";
         if(value==test){
            program::pump_function=program::two_temperature;
            return true;
         }
         test="double-pulse-two-temperature";
         if(value==test){
            program::pump_function=program::double_pump_two_temperature;
            return true;
         }
         test="double-pulse-square";
         if(value==test){
            program::pump_function=program::double_pump_square;
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
            std::cerr << "\t\"square\"" << std::endl;
            std::cerr << "\t\"double-pulse-square\"" << std::endl;
            std::cerr << "\t\"two-temperature\"" << std::endl;
            std::cerr << "\t\"double-pulse-two-temperature\"" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }
      //--------------------------------------------------------------------
      test="laser-pulse-time";
      if(word==test){
         double pt=atof(value.c_str());
         vin::check_for_valid_value(pt, word, line, prefix, unit, "time", 1.0e-18, 1.0,"input","1 attosecond - 1 s");
         program::pump_time=pt;
         return true;
      }
      //--------------------------------------------------------------------
      test="laser-pulse-power";
      if(word==test){
         double pp=atof(value.c_str());
         vin::check_for_valid_value(pp, word, line, prefix, unit, "none", 0.0, 1.0e40,"input","0.0 - 1.0E40");
         program::pump_power=pp;
         return true;
      }
      //--------------------------------------------------------------------
      test="second-laser-pulse-time";
      if(word==test){
         double pt=atof(value.c_str());
         vin::check_for_valid_value(pt, word, line, prefix, unit, "time", 1.0e-18, 1.0,"input","1 attosecond - 1 s");
         program::double_pump_time=pt;
         return true;
      }
      //--------------------------------------------------------------------
      test="second-laser-pulse-power";
      if(word==test){
         double pp=atof(value.c_str());
         vin::check_for_valid_value(pp, word, line, prefix, unit, "none", 0.0, 1.0e40,"input","0.0 - 1.0E40");
         program::double_pump_power=pp;
         return true;
      }
      //--------------------------------------------------------------------
      test="second-laser-pulse-maximum-temperature";
      if(word==test){
         double T=atof(value.c_str());
         vin::check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 1.0e6,"input","0.0 - 1,000,000 K");
         program::double_pump_Tmax=T;
         return true;
      }
      //--------------------------------------------------------------------
      test="second-laser-pulse-delay-time";
      if(word==test){
         double pd=atof(value.c_str());
         vin::check_for_valid_value(pd, word, line, prefix, unit, "time", 0.0, 1.0,"input","0 - 1 s");
         program::double_pump_delay=pd;
         return true;
      }
      //--------------------------------------------------------------------
      test="two-temperature-heat-sink-coupling";
      if(word==test){
         double hscc=atof(value.c_str());
         vin::check_for_valid_value(hscc, word, line, prefix, unit, "none", 0.0, 1.0e40,"input","0.0 - 1.0E40");
         program::HeatSinkCouplingConstant=hscc;
         return true;
      }
      //--------------------------------------------------------------------
      test="two-temperature-electron-heat-capacity";
      if(word==test){
         double hscc=atof(value.c_str());
         vin::check_for_valid_value(hscc, word, line, prefix, unit, "none", 0.0, 1.0e40,"input","0.0 - 1.0E40");
         program::TTCe=hscc;
         return true;
      }
      //--------------------------------------------------------------------
      test="two-temperature-phonon-heat-capacity";
      if(word==test){
         double hscc=atof(value.c_str());
         vin::check_for_valid_value(hscc, word, line, prefix, unit, "none", 0.0, 1.0e40,"input","0.0 - 1.0E40");
         program::TTCl=hscc;
         return true;
      }
      //--------------------------------------------------------------------
      test="two-temperature-electron-phonon-coupling";
      if(word==test){
         double hscc=atof(value.c_str());
         vin::check_for_valid_value(hscc, word, line, prefix, unit, "none", 0.0, 1.0e40,"input","0.0 - 1.0E40");
         program::TTG=hscc;
         return true;
      }
      //--------------------------------------------------------------------
      test="applied-field-strength";
      if(word==test){
         double H=atof(value.c_str());
         vin::check_for_valid_value(H, word, line, prefix, unit, "field", -1.e4, 1.0e4,"input","+/- 10,000 T");
         sim::H_applied=H;
         return true;
      }
      //--------------------------------------------------------------------
      test="minimum-applied-field-strength";
      if(word==test){
         double H=atof(value.c_str());
         vin::check_for_valid_value(H, word, line, prefix, unit, "field", 0.0, 1.0e3,"input","0 - 1,000 T");
         sim::Hmin=H;
         return true;
      }
      //--------------------------------------------------------------------
      test="maximum-applied-field-strength";
      if(word==test){
         double H=atof(value.c_str());
         vin::check_for_valid_value(H, word, line, prefix, unit, "field", 0.0, 1.0e3,"input","0 - 1,000 T");
         sim::Hmax=H;
         return true;
      }
      //--------------------------------------------------------------------
      test="equilibration-applied-field-strength";
      if(word==test){
         double H=atof(value.c_str());
         vin::check_for_valid_value(H, word, line, prefix, unit, "field", 0.0, 1.0e3,"input","0 - 1,000 T");
         sim::Heq=H;
         return true;
      }
      //--------------------------------------------------------------------
      test="applied-field-strength-increment";
      if(word==test){
         double H=atof(value.c_str());
         vin::check_for_valid_value(H, word, line, prefix, unit, "field", 1.0e-6, 1.0e3,"input","1 uT - 1,000 T");
         sim::Hinc=H;
         return true;
      }
      //--------------------------------------------------------------------
      test="applied-field-angle-theta";
      if(word==test){
         double angle=atof(value.c_str());
         vin::check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
         sim::applied_field_angle_theta=angle;
         sim::applied_field_set_by_angle=true;
         return true;
      }
      //--------------------------------------------------------------------
      test="applied-field-angle-phi";
      if(word==test){
         double angle=atof(value.c_str());
         vin::check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
         sim::applied_field_angle_phi=angle;
         sim::applied_field_set_by_angle=true;
         return true;
      }
      //--------------------------------------------------------------------
      test="applied-field-unit-vector";
      if(word==test){
         std::vector<double> u(3);
         u=vin::doubles_from_string(value);
         vin::check_for_valid_unit_vector(u, word, line, prefix, "input");
         sim::H_vec[0]=u.at(0);
         sim::H_vec[1]=u.at(1);
         sim::H_vec[2]=u.at(2);
         sim::applied_field_set_by_angle=false;
         return true;
      }
      //--------------------------------------------------------------------
      test="demagnetisation-factor";
      if(word==test){
         std::vector<double> u(3);
         u=vin::doubles_from_string(value);
         vin::check_for_valid_three_vector(u, word, line, prefix, "input");
         // Extra check for demagnetisation-factor Nx+Ny+Nz=1
         double sum=u.at(0)+u.at(1)+u.at(2);
         if(fabs(1.0-sum)>1.e-4){
            std::ostringstream message;
            message << "Error: sum of all elements of variable " << prefix << word << " on line " << line << " of input file must equal 1.";
            err::fatal(message.str());
         }
         sim::demag_factor[0]=u.at(0);
         sim::demag_factor[1]=u.at(1);
         sim::demag_factor[2]=u.at(2);
         sim::ext_demag=true;
         // force calculation of system magnetization
         stats::calculate_system_magnetization=true;
         return true;
      }
      //--------------------------------------------------------------------
      test="fmr-field-strength";
      if(word==test){
         double H=atof(value.c_str());
         vin::check_for_valid_value(H, word, line, prefix, unit, "field", -1.e4, 1.0e4,"input","+/- 10,000 T");
         sim::fmr_field_strength=H;
         return true;
      }
      //--------------------------------------------------------------------
      test="fmr-field-frequency";
      if(word==test){
         double w = atof(value.c_str());
         vin::check_for_valid_value(w, word, line, prefix, unit, "frequency", 0.0, 1.0e14,"input","0 - 10,000 GHz");
         sim::fmr_field_frequency = w;
         return true;
      }
      //--------------------------------------------------------------------
      test="fmr-field-unit-vector";
      if(word==test){
         std::vector<double> u(3);
         u=vin::doubles_from_string(value);
         vin::check_for_valid_unit_vector(u, word, line, prefix, "input");
         sim::fmr_field_unit_vector = u;
         return true;
      }
      //--------------------------------------------------------------------
      test="constraint-rotation-update";
      if(word==test){
         montecarlo::constraint_rotation=true;
         return true;
      }
      //--------------------------------------------------------------------
      test="constraint-angle-theta";
      if(word==test){
         double angle=atof(value.c_str());
         vin::check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
         montecarlo::constraint_theta=angle;
         return true;
      }
      //--------------------------------------------------------------------
      test="constraint-angle-theta-minimum";
      if(word==test){
         double angle=atof(value.c_str());
         vin::check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
         montecarlo::constraint_theta_min=angle;
         return true;
      }
      //--------------------------------------------------------------------
      test="constraint-angle-theta-maximum";
      if(word==test){
         double angle=atof(value.c_str());
         vin::check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
         montecarlo::constraint_theta_max=angle;
         return true;
      }
      //--------------------------------------------------------------------
      test="constraint-angle-theta-increment";
      if(word==test){
         double angle=atof(value.c_str());
         vin::check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
         montecarlo::constraint_theta_delta=angle;
         return true;
      }
      //--------------------------------------------------------------------
      test="constraint-angle-phi";
      if(word==test){
         double angle=atof(value.c_str());
         vin::check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
         montecarlo::constraint_phi=angle;
         return true;
      }
      //--------------------------------------------------------------------
      test="constraint-angle-phi-minimum";
      if(word==test){
         double angle=atof(value.c_str());
         vin::check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
         montecarlo::constraint_phi_min=angle;
         return true;
      }
      //--------------------------------------------------------------------
      test="constraint-angle-phi-maximum";
      if(word==test){
         double angle=atof(value.c_str());
         vin::check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
         montecarlo::constraint_phi_max=angle;
         return true;
      }
      //--------------------------------------------------------------------
      test="constraint-angle-phi-increment";
      if(word==test){
         double angle=atof(value.c_str());
         vin::check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
         montecarlo::constraint_phi_delta=angle;
         return true;
      }
      //--------------------------------------------------------------------
      test="save-checkpoint";
      if(word==test){
         test="end";
         if(value==test){
            sim::save_checkpoint_flag=true; // Save checkpoint
            sim::save_checkpoint_continuous_flag=false; // do not save checkpoints during simulation
            return true;
         }
         test="continuous";
         if(value==test){
            sim::save_checkpoint_flag=true; // Save checkpoint
            sim::save_checkpoint_continuous_flag=true; // save checkpoints during simulation
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
            std::cerr << "\t\"end\"" << std::endl;
            std::cerr << "\t\"continuous\"" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }
      //--------------------------------------------------------------------
      test="save-checkpoint-rate";
      if(word==test){
         int scr=atoi(value.c_str());
         vin::check_for_valid_int(scr, word, line, prefix, 1, 2000000000,"input","1 - 2,000,000,000");
         sim::save_checkpoint_rate=scr;
         return true;
      }
      //--------------------------------------------------------------------
      test="load-checkpoint";
      if(word==test){
         test="restart";
         if(value==test){
            sim::load_checkpoint_flag=true; // Load spin configurations
            sim::load_checkpoint_continue_flag=false; // Restart simulation with checkpoint configuration
            return true;
         }
         test="continue";
         if(value==test){
            sim::load_checkpoint_flag=true; // Load spin configurations
            sim::load_checkpoint_continue_flag=true; // Continue simulation from saved time with checkpoint configuration
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
            std::cerr << "\t\"restart\"" << std::endl;
            std::cerr << "\t\"continue\"" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }
      //--------------------------------------------------------------------
      test="load-checkpoint-if-exists";
      if(word==test){
         // determine checkpoint file name
         std::stringstream chkfilenamess;
         chkfilenamess << "vampire" << vmpi::my_rank << ".chk";
         std::string chkfilename = chkfilenamess.str();

         // open checkpoint file
         std::ifstream chkfile(chkfilename.c_str());
         bool exists=chkfile.good();
         chkfile.close();
         if(exists){
            test="restart";
            if(value==test){
               sim::load_checkpoint_flag=true; // Load spin configurations
               sim::load_checkpoint_continue_flag=false; // Restart simulation with checkpoint configuration
               return true;
            }
            test="continue";
            if(value==test){
               sim::load_checkpoint_flag=true; // Load spin configurations
               sim::load_checkpoint_continue_flag=true; // Continue simulation from saved time with checkpoint configuration
               return true;
            }
            else{
               terminaltextcolor(RED);
               std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
               std::cerr << "\t\"restart\"" << std::endl;
               std::cerr << "\t\"continue\"" << std::endl;
               terminaltextcolor(WHITE);
               err::vexit();
            }
         }
         else{
            return true; // No checkpoint file
         }
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
      std::string test="domain-wall-second-magnetisation-vector";
      if(word==test){
         std::vector<double> u(3);
         u=vin::doubles_from_string(value);
         vin::check_for_valid_unit_vector(u, word, line, prefix, "input");
         std::cout << program::domain_wall_second_vector_x.size() << "\t" << super_index << "\t" << u[0] << '\t' << u[1] << '\t' << u[2] <<std::endl;
         program::domain_wall_second_vector_x[super_index] = u[0];
         program::domain_wall_second_vector_y[super_index] = u[1];
         program::domain_wall_second_vector_z[super_index] = u[2];
         return true;
      }
      //------------------------------------------------------------
      test  = "spin-transfer-relaxation-torque";
      std::string test2 = "slonczewski-adiabatic-spin-torque";
      std::string test3 = "spin-transfer-torque";
      std::string test4 = "antidamping-torque";
      // aj parameter for material in slonczewski torque calculation
      if( word==test || word==test2 || word==test3 || word==test4){
         double aj=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(aj, word, line, prefix, unit, "field", -1.0e-2, 1.0e2,"input","-100 - 100T");
         sim::internal::mp[super_index].stt_rj.set(aj);
         sim::internal::enable_spin_torque_fields = true;
         return true;
      }
      //------------------------------------------------------------
      test2 = "spin-transfer-precession-torque";
      test  = "slonczewski-non-adiabatic-spin-torque";
      test3 = "field-like-torque";
      test4 = "slonczewski-precession-spin-torque";
      // bj parameter for material in slonczewski torque calculation
      if( word==test || word==test2 || word==test3 ){
         double bj=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(bj, word, line, prefix, unit, "field", -1.0e-2, 1.0e2,"input","-100 - 100T");
         sim::internal::mp[super_index].stt_pj.set(bj);
         sim::internal::enable_spin_torque_fields = true;
         return true;
      }
      //------------------------------------------------------------
      test = "spin-transfer-torque-asymmetry";
      // damping-like parameter for material in spin orbit torque calculation
      if( word==test ){
         double sttasm = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(sttasm, word, line, prefix, unit, "", 0.0, 1.0e2,"input","0 - 100");
         sim::internal::mp[super_index].stt_asm.set(sttasm);
         return true;
      }
      //------------------------------------------------------------
      test="spin-transfer-torque-polarization-unit-vector";
      if(word==test){
         std::vector<double> u(3);
         u=vin::doubles_from_string(value);
         // Test for valid range
         vin::check_for_valid_unit_vector(u, word, line, prefix, "input");
         // save sanitized unit vector
         sim::internal::mp[super_index].stt_pv.x = u[0];
         sim::internal::mp[super_index].stt_pv.y = u[1];
         sim::internal::mp[super_index].stt_pv.z = u[2];
         sim::internal::enable_local_stt_polarizers = true;
         return true;
      }
      //------------------------------------------------------------
      test="spin-transfer-torque-polarization-material";
      if(word==test){
         int mid = vin::str_to_int(value); // convert string to int
         // Test for valid range
         vin::check_for_valid_int(mid, word, line, prefix, 1, mp::max_materials,"material","1 - 100");
         sim::internal::mp[super_index].stt_pm.set(mid-1); // save material ID as m-1
         zlog << zTs() << "Setting STT polarizer for material " << super_index+1 << " to material " << sim::internal::mp[super_index].stt_pm.get()+1 << std::endl;
         sim::internal::enable_local_stt_polarizers = true;
         return true;
      }
      //------------------------------------------------------------
      // field-like parameter for material in spin orbit torque calculation
      test = "spin-orbit-relaxation-torque";
      test2 = "spin-orbit-anti-damping-torque";
      if( word==test || word==test2 ){
         double aj = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(aj, word, line, prefix, unit, "field", -1.0e2, 1.0e2,"input","-100 - 100T");
         sim::internal::mp[super_index].sot_rj.set(aj);
         sim::internal::enable_spin_torque_fields = true;
         return true;
      }
      //------------------------------------------------------------
      test = "spin-orbit-precession-torque";
      test2 = "spin-orbit-torque";
      test3 = "spin-orbit-field-like-torque";
      // damping-like parameter for material in spin orbit torque calculation
      if( word==test || word==test2 || word==test3 ){
         double bj = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(bj, word, line, prefix, unit, "field", -1.0e2, 1.0e2,"input","-100 - 100T");
         sim::internal::mp[super_index].sot_pj.set(bj);
         sim::internal::enable_spin_torque_fields = true;
         return true;
      }
      //------------------------------------------------------------
      test = "spin-orbit-torque-asymmetry";
      // damping-like parameter for material in spin orbit torque calculation
      if( word==test ){
         double sotasm = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(sotasm, word, line, prefix, unit, "", 0.0, 1.0e2,"input","0 - 100");
         sim::internal::mp[super_index].sot_asm.set(sotasm);
         return true;
      }
      //------------------------------------------------------------
      test="spin-orbit-torque-polarization-unit-vector";
      if(word==test){
         std::vector<double> u(3);
         u=vin::doubles_from_string(value);
         // Test for valid range
         vin::check_for_valid_unit_vector(u, word, line, prefix, "input");
         // save sanitized unit vector
         sim::internal::mp[super_index].sot_pv.x = u[0];
         sim::internal::mp[super_index].sot_pv.y = u[1];
         sim::internal::mp[super_index].sot_pv.z = u[2];
         sim::internal::enable_local_sot_polarizers = true;
         return true;
      }
      //------------------------------------------------------------
      test="spin-orbit-torque-polarization-material";
      if(word==test){
         int mid = vin::str_to_int(value); // convert string to int
         // Test for valid range
         vin::check_for_valid_int(mid, word, line, prefix, 1, mp::max_materials,"material","1 - 100");
         sim::internal::mp[super_index].sot_pm.set(mid-1); // save material ID as m-1
         zlog << zTs() << "Setting SOT polarizer for material " << super_index+1 << " to material " << sim::internal::mp[super_index].sot_pm.get()+1 << std::endl;
         sim::internal::enable_local_sot_polarizers = true;
         return true;
      }
      //------------------------------------------------------------
      test  = "voltage-controlled-magnetic-anisotropy-coefficient";
      test2 = "vcma-coefficient";
      // voltage dependent anisotropy (Joules / Volt == Coulombs) for vcma
      if( word == test || word == test2 ){
         double vcmak = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(vcmak, word, line, prefix, unit, "", -1.0e-15, 1.0e-15,"input","0 - 1e-21");
         sim::internal::mp[super_index].vcmak.set(vcmak);
         sim::internal::enable_vcma_fields = true; // enable vcma fields
         return true;
      }
      //--------------------------------------------------------------------
      test="temperature";
      if(word==test){
         double T=atof(value.c_str());
         // Test for valid range
         if((T>=0.0) && (T<1.0E5)){
            mp::read_material[super_index].temperature=T;
            // set local temperature flag
            sim::local_temperature=true;
            return true;
         }
         else{
            std::ostringstream message;
            message << "Error on line " << line << " of material file - material[" << super_index+1 << "]:" << word << " is outside of valid range 0.0 - 1.0E5";
            err::fatal(message.str());
         }
      }
      //--------------------------------------------------------------------
      test="minimum-temperature";
      if(word==test){
         double T=atof(value.c_str());
         vin::check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 1.0e6,"input","0.0 - 1,000,000 K");
         mp::read_material[super_index].minimum_temperature = T;
         sim::local_temperature=true; // set local temperature flag
         return true;
      }
      //--------------------------------------------------------------------
      test="maximum-temperature";
      if(word==test){
         double T=atof(value.c_str());
         vin::check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 1.0e6,"input","0.0 - 1,000,000 K");
         mp::read_material[super_index].maximum_temperature = T;
         sim::local_temperature=true; // set local temperature flag
         return true;
      }
      //--------------------------------------------------------------------
      test="use-phonon-temperature";
      /*
      logical use-phonon-temperature
          This flag enables specific materials to couple to the phonon temperature
          of the system for simulations using the two temperature model. The default
          is for all materials to use the electron temperature. Valid values are true,
          false or (blank) [same as true].
      */
      if(word==test){
         // Test for sane input
         //bool sanitised_bool=vin::check_for_valid_bool(value, word, line, prefix,"material");

         // set flag
         mp::read_material[super_index].couple_to_phonon_temperature=true; //sanitised_bool;

         // enable local temperature flag
         sim::local_temperature=true;
         return true;
      }
      //--------------------------------------------------------------------
      test="applied-field-strength";
      if(word==test){
         double H=atof(value.c_str());
         // test for unit
         std::string unit_type="field";
         // if no unit given, assume internal
         if(unit.size() != 0){
            units::convert(unit,H,unit_type);
         }
         std::string str="field";
         if(unit_type==str){
            // Test for valid range
            if((H>=-1.0E5) && (H<1.0E5)){
               mp::read_material[super_index].applied_field_strength=H;
               // set local applied field flag
               sim::local_applied_field=true;
               return true;
            }
            else{
               std::ostringstream message;
               message << "Error on line " << line << " of material file - material[" << super_index+1 << "]:" << word << " is outside of valid range -1.0e5 - +1.0e5";
               err::fatal(message.str());
            }
         }
         else{
            std::ostringstream message;
            message << "Error on line " << line << " of material file - unit type \'" << unit_type << "\' is invalid for parameter material[" << super_index+1 << "]:" << word;
            err::fatal(message.str());
         }
      }
      //--------------------------------------------------------------------
      test="applied-field-unit-vector";
      if(word==test){
         // temporary storage container
         std::vector<double> u(3);

         // read values from string
         u=vin::doubles_from_string(value);

         // check size
         if(u.size()!=3){
            std::ostringstream message;
            message << "Error in input file - material[" << super_index+1 << "]:" << word << " must have three values.";
            err::fatal(message.str());
         }

         // Normalise
         double ULength=sqrt(u.at(0)*u.at(0)+u.at(1)*u.at(1)+u.at(2)*u.at(2));

         // Check for correct length unit vector
         if(ULength < 1.0e-9){
            std::ostringstream message;
            message << "Error in input file - material[" << super_index+1 << "]:" << word << " must be normalisable (possibly all zero).";
            err::fatal(message.str());
         }
         u.at(0)/=ULength;
         u.at(1)/=ULength;
         u.at(2)/=ULength;

         // Copy applied field direction to material
         mp::read_material[super_index].applied_field_unit_vector=u;

         // set local applied field flag
         sim::local_applied_field=true;

         return true;

      }
      //--------------------------------------------------------------------
      test="fmr-field-strength";
      if(word==test){
         double H = atof(value.c_str());
         vin::check_for_valid_value(H, word, line, prefix, unit, "field", 0.0, 1.0e5,"material","0 - 10,000 T");
         mp::read_material[super_index].fmr_field_strength=H;
         // set local fmr flag
         sim::local_fmr_field=true;
         return true;
      }
      //--------------------------------------------------------------------
      test="fmr-field-frequency";
      if( word == test ){
         double f = atof(value.c_str());
         vin::check_for_valid_value(f, word, line, prefix, unit, "frequency", 0.0, 1.0e14,"material","0 - 10,000 GHz");
         mp::read_material[super_index].fmr_field_frequency = f;
         sim::local_fmr_field = true;
         return true;
      }
      //--------------------------------------------------------------------
      test="fmr-field-unit-vector";
      if(word==test){
         // temporary storage container
         std::vector<double> u(3);

         // read values from string
         u=vin::doubles_from_string(value);

         // check size
         if(u.size()!=3){
            std::ostringstream message;
            message << "Error on line " << line << " of material file - material[" << super_index+1 << "]:" << word << " must have three values.";
            err::fatal(message.str());
         }

         // Normalise
         double ULength=sqrt(u.at(0)*u.at(0)+u.at(1)*u.at(1)+u.at(2)*u.at(2));

         // Check for correct length unit vector
         if(ULength < 1.0e-9){
            std::ostringstream message;
            message << "Error on line " << line << " of material file - material[" << super_index+1 << "]:" << word << " must be normalisable (possibly all zero).";
            err::fatal(message.str());
         }
         u.at(0)/=ULength;
         u.at(1)/=ULength;
         u.at(2)/=ULength;

         // Copy field direction to material
         mp::read_material[super_index].fmr_field_unit_vector=u;

         // set local applied field flag
         sim::local_fmr_field=true;

         return true;

      }
      //--------------------------------------------------------------------
      // keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of namespace sim
