//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2009-2012. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
///
/// @file
/// @brief Contains the Temperature Pulse program
///
/// @details Performs a time series with variable (pulsed) temperature
///
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 2.0
/// @date    08/11/2012
///=====================================================================================
///

// Standard Libraries
#include <cmath>
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "constants.hpp"
#include "errors.hpp"
#include "gpu.hpp"
#include "material.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

/// function forward declarations
double two_temperature_function(double ftime);
double double_pump_two_temperature_function(double ftime);
double square_temperature_function(double ftime);
double double_pump_square_temperature_function(double ftime);
double multilayer_two_temperature_function(double ftime);

/// Selects appropriate temperature pulse function according to user input
double temperature_pulse_function(double function_time){

   switch(program::pump_function){

      case(program::square):
         return square_temperature_function(function_time);
         break;

      case(program::two_temperature):
         return two_temperature_function(function_time);
         break;

      case(program::double_pump_two_temperature):
         return double_pump_two_temperature_function(function_time);
         break;

      case(program::double_pump_square):
         return double_pump_square_temperature_function(function_time);
         break;

      case(program::multilayer_two_temperature):
         return multilayer_two_temperature_function(function_time);
         break;

      default:
         return 0.0;
   }
}

///-----------------------------------------------------------------------------------------
///   Calculates temperature for single Gaussian pulse using the two-temperature model
///
///   Uses the free electron approximation in calculating the electron temperature
///   where C_e = gamma_e T_e (ref Unai Atxitia PhD thesis pp 53).
///   gamma_e ~ 3e3 J/m^-3/K^-2
///
///-----------------------------------------------------------------------------------------
///   Pump power estimation and units
///
///   Improved Two-Temperature Model and Its Application in Ultrashort Laser Heating of Metal Films
///   Lan Jiang and Hai-Lung Tsai J. Heat Transfer 127(10), 1167-1173 (Jun 07, 2005) (7 pages)
///
///   General TT model has a source term S(z,t) = P I(z,t) where P is the laser power and I is
///   the laser intensity. Assuming a penetration depth delta gives
///
///   I(z,t) = 2/(tp delta sqrt(pi/ln 2)) (1-R) exp[-(4 ln 2)(t/tp)**2] exp[-z/delta] W/m^3
///
///   where z is the height from the surface, tp is the laser pulse time, R is the reflectivity
///   coefficient, and t is the time. Assuming all energy from the laser is absorbed in the
///   sample (R = 1, z >> delta) and uniform heating (infinite thermal conductivity) this leads
///   to a simplified expression for the laser intensity
///
///   I(t) = 2/(tp delta sqrt(pi/ln 2)) exp[-(4 ln 2)(t/tp)**2] W/m^3
///
///   Assuming a nominal penetration depth of 10 nm and conversion from mJ/cm^2 -> J/m^2,
///   I(t) is now in units of s^1 m^-1 and laser power P is in units of mJ/cm^2
///
///-----------------------------------------------------------------------------------------
double two_temperature_function(double ftime){

   const double i_pump_time = 1.0/program::pump_time;
   const double reduced_time = (ftime-3.*program::pump_time)*i_pump_time;
   const double four_ln_2 = 2.77258872224; // 4 ln 2
   // 2/(delta sqrt(pi/ln 2))*10, delta = 10 nm, mJ/cm^2 -> J/m^2 (factor 10)
   const double two_delta_sqrt_pi_ln_2 = 939437278.7;
   const double pump=program::pump_power*two_delta_sqrt_pi_ln_2*
                     exp(-four_ln_2*reduced_time*reduced_time)*i_pump_time;

   const double Te = program::TTTe;
   const double Tp = program::TTTp;
   const double G  = program::TTG;
   const double Ce = program::TTCe;
   const double Cl = program::TTCl;
   const double dt = mp::dt_SI;

   // integrate two temperature model (floor in free elecron approximation (c prop to T) for low temperatures)
   if(Te>1.0) program::TTTe = (-G*(Te-Tp)+pump)*dt/(Ce*Te) + Te;
   else program::TTTe =       (-G*(Te-Tp)+pump)*dt/Ce + Te;
   program::TTTp =            ( G*(Te-Tp)     )*dt/Cl + Tp - (Tp-sim::Teq)*program::HeatSinkCouplingConstant*dt;

   // Optionally set material specific temperatures
   if(sim::local_temperature==true){
      for(unsigned int mat=0;mat<mp::material.size();mat++){
         if(mp::material[mat].couple_to_phonon_temperature==true) mp::material[mat].temperature=program::TTTp;
         else mp::material[mat].temperature=program::TTTe;
      }
   }

   return program::TTTe;

}

///-----------------------------------------------------------------------------------------
/// Multilayer two-temperature model (M Strungaru)
///
///   Until sim:two-temperature-cutoff-time the film is heated by the same Gaussian pulse as
///   two_temperature_function, with the film Ce, Cl and G and no heat-sink coupling, and the
///   electron temperature is returned. From the cutoff time the two-temperature model uses
///   the substrate parameters (two-temperature-substrate-*) with heat-sink coupling and the
///   phonon temperature is returned. If sim:two-temperature-damping is set, the Gilbert
///   damping of every material is changed to it at the cutoff time.
///
///-----------------------------------------------------------------------------------------
static bool multilayer_damping_switched = false;

static void set_damping_of_all_materials(const double alpha){
   for(unsigned int mat=0;mat<mp::material.size();mat++){
      mp::material[mat].alpha = alpha;
      mp::material[mat].one_oneplusalpha_sq   = -mp::material[mat].gamma_rel/(1.0+alpha*alpha);
      mp::material[mat].alpha_oneplusalpha_sq =  alpha*mp::material[mat].one_oneplusalpha_sq;
      mp::material[mat].H_th_sigma            = sqrt(2.0*alpha*constants::kB / (mp::material[mat].mu_s_SI*mp::material[mat].gamma_rel*mp::dt));
   }
}

double multilayer_two_temperature_function(double ftime){

   // material parameters are copied to the GPU at initialisation, so the damping switch would be lost
   if(gpu::acceleration && program::TTdamp >= 0.0){
      err::fatal("Error: sim:two-temperature-damping with sim:laser-pulse-temporal-profile = multilayer-two-temperature is not implemented on the GPU.");
   }

   const double i_pump_time = 1.0/program::pump_time;
   const double reduced_time = (ftime-3.*program::pump_time)*i_pump_time;
   const double four_ln_2 = 2.77258872224; // 4 ln 2
   // 2/(delta sqrt(pi/ln 2))*10, delta = 10 nm, mJ/cm^2 -> J/m^2 (factor 10)
   const double two_delta_sqrt_pi_ln_2 = 939437278.7;
   const double pump=program::pump_power*two_delta_sqrt_pi_ln_2*
                     exp(-four_ln_2*reduced_time*reduced_time)*i_pump_time;

   const double Te = program::TTTe;
   const double Tp = program::TTTp;
   const double dt = mp::dt_SI;

   const bool after_cutoff = ftime >= program::TTcutofftime;

   if(after_cutoff){

      // change the damping once at the cutoff time
      if(!multilayer_damping_switched && program::TTdamp >= 0.0) set_damping_of_all_materials(program::TTdamp);
      multilayer_damping_switched = true;

      // substrate parameters with heat-sink coupling
      const double G  = program::TTG_sub;
      const double Ce = program::TTCe_sub;
      const double Cl = program::TTCl_sub;
      if(Te>1.0) program::TTTe = (-G*(Te-Tp)+pump)*dt/(Ce*Te) + Te;
      else program::TTTe =       (-G*(Te-Tp)+pump)*dt/Ce + Te;
      program::TTTp =            ( G*(Te-Tp)     )*dt/Cl + Tp - (Tp-sim::Teq)*program::HeatSinkCouplingConstant*dt;

   }
   else{

      // film parameters without heat-sink coupling
      const double G  = program::TTG;
      const double Ce = program::TTCe;
      const double Cl = program::TTCl;
      if(Te>1.0) program::TTTe = (-G*(Te-Tp)+pump)*dt/(Ce*Te) + Te;
      else program::TTTe =       (-G*(Te-Tp)+pump)*dt/Ce + Te;
      program::TTTp =            ( G*(Te-Tp)     )*dt/Cl + Tp;

   }

   // Optionally set material specific temperatures
   if(sim::local_temperature==true){
      for(unsigned int mat=0;mat<mp::material.size();mat++){
         if(mp::material[mat].couple_to_phonon_temperature==true) mp::material[mat].temperature=program::TTTp;
         else mp::material[mat].temperature=program::TTTe;
      }
   }

   return after_cutoff ? program::TTTp : program::TTTe;

}

/// Calculates temperature for double Gaussian pulses using the two-temperature model
double double_pump_two_temperature_function(double ftime){

   const double four_ln_2 = 2.77258872224; // 4 ln 2
   // 2/(delta sqrt(pi/ln 2))*10, delta = 10 nm, mJ/cm^2 -> J/m^2 (factor 10)
   const double two_delta_sqrt_pi_ln_2 = 939437278.7;

   const double i_pump_time1 = 1.0/program::pump_time;
   const double reduced_time1 = (ftime-3.*program::pump_time)*i_pump_time1;
   const double pump1=program::pump_power*two_delta_sqrt_pi_ln_2*
                     exp(-four_ln_2*reduced_time1*reduced_time1)*i_pump_time1;

   const double i_pump_time2 = 1.0/program::double_pump_time;
   const double reduced_time2 = (ftime-program::double_pump_delay-3.*program::double_pump_time)*i_pump_time2;
   const double pump2=program::double_pump_power*two_delta_sqrt_pi_ln_2*
                     exp(-four_ln_2*reduced_time2*reduced_time2)*i_pump_time2;

      const double Te = program::TTTe;
      const double Tp = program::TTTp;
      const double G  = program::TTG;
      const double Ce = program::TTCe;
      const double Cl = program::TTCl;
      const double dt = mp::dt_SI;

      // integrate two temperature model (floor in free elecron approximation (c prop to T) for low temperatures)
      if(Te>1.0) program::TTTe = (-G*(Te-Tp)+pump1+pump2)*dt/(Ce*Te) + Te;
      else program::TTTe =       (-G*(Te-Tp)+pump1+pump2)*dt/Ce + Te;
      program::TTTp =            ( G*(Te-Tp)            )*dt/Cl + Tp - (Tp-sim::Teq)*program::HeatSinkCouplingConstant*dt;

      // Optionally set material specific temperatures
      if(sim::local_temperature==true){
         for(unsigned int mat=0;mat<mp::material.size();mat++){
            if(mp::material[mat].couple_to_phonon_temperature==true) mp::material[mat].temperature=program::TTTp;
            else mp::material[mat].temperature=program::TTTe;
         }
      }

      return program::TTTe;
}

/// Calculates temperature for single square pulse
double square_temperature_function(double ftime){

   //
   ///
   ///
   ///            Tmax   |---------------|
   ///                   |               |
   ///                   |               |
   /// ------------------|               |--------------- Teq
   ///                    <- pump time ->
   //

   if(ftime <= program::pump_time) return sim::Tmax;
   else return sim::Teq;

}

/// Calculates temperature for double square pulse
double double_pump_square_temperature_function(double ftime){

   ///
   ///          0               t1                t2              t3
   ///
   ///   Tmax   |---------------|
   ///          |               |                 |---------------| double_pump_Tmax
   ///          |               |                 |               |
   /// ---------|  -  -  -  -   |-----------------|   -   -   -   |------------ Teq
   ///
   ///           <- pump time -> <- dpump_delay -> <- pump time2 ->
   ///

   const double t1 = program::pump_time;
   const double t2 = t1+program::double_pump_delay;
   const double t3 = t2+program::double_pump_time;

   if(ftime <= t1) return sim::Tmax;
   else if( (ftime > t2) && (ftime <=t3) ) return program::double_pump_Tmax;
   else return sim::Teq;

}

namespace program{
/// @brief Function to calculate a time series with two temperature model heating/cooling profile
///
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.1
/// @date    06/04/2011
///
/// @internal
///   Created:    10/03/2010
///   Revision:   06/04/2011
///=====================================================================================
///
void temperature_pulse(){

   // check calling of routine if error checking is activated
   if(err::check==true){std::cout << "program::two_temperature_pulse has been called" << std::endl;}

   // Set equilibration temperature and field
   sim::temperature=sim::Teq;

   // Initialise electron and phonon temperature
   program::TTTe=sim::temperature;
   program::TTTp=sim::temperature;

   // If local temperature is set then also initalise local temperatures
   if(sim::local_temperature==true){
      for(unsigned int mat=0;mat<mp::material.size();mat++){
         if(mp::material[mat].couple_to_phonon_temperature==true) mp::material[mat].temperature=program::TTTp;
         else mp::material[mat].temperature=program::TTTe;
      }
   }

   // Equilibrate system
   while(sim::time<sim::equilibration_time){

      sim::integrate(sim::partial_time);

      // Calculate magnetisation statistics
      stats::update();

      // Output data
      vout::data();
   }

   //loop sim::runs times
   for(int r=0; r<sim::runs;r++){

   // record starting time after equiibration/last pulse
   uint64_t start_time = sim::time;

   // Simulate temperature pulse
   while(sim::time<sim::total_time+start_time){

      // loop over partial_time to update temperature every time
      for(uint64_t tt=0; tt < sim::partial_time; tt++){

         // Calculate time from pulse
         double time_from_start=mp::dt_SI*double(sim::time-start_time);

         // Calculate temperature
         sim::temperature=temperature_pulse_function(time_from_start);

         // Integrate system
         sim::integrate(1);

      }

      // Calculate magnetisation statistics
      stats::update();

      // Output data
      vout::data();

   }

   }

} // end of temperature pulse

} // end of namespace program
