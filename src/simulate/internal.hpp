#ifndef SIM_INTERNAL_H_
#define SIM_INTERNAL_H_
//-----------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

//---------------------------------------------------------------------
// Defines shared internal data structures and functions for the
// simulation methods implementation. These functions should
// not be accessed outside of the simulate module.
//---------------------------------------------------------------------

namespace sim{
   namespace internal{

      //-----------------------------------------------------------------------------
      // Internal shared variables used for the simulation
      //-----------------------------------------------------------------------------
      // simple initialised class for set variables
      class set_double_t{

      private:
         double value; // value
         bool setf; // flag specifiying variable has been set

      public:
         // class functions
         // constructor
         set_double_t() : value(0.0), setf(false) { }

         // setting function
         void set(double in_value){
            value = in_value;
            setf = true;
         };

         // get value function
         double get(){ return value; };
         // check if variable is set
         bool is_set(){ return setf; };

      };

      struct mp_t{
         set_double_t stt_rj; // spin tranfer relaxation torque
         set_double_t stt_pj; // spin transfer precession torque
         set_double_t sot_rj; // spin orbit relaxation torque
         set_double_t sot_pj; // spin orbit precession torque
      };

      extern bool enable_spin_torque_fields; // flag to enable spin torque fields

      extern std::vector<sim::internal::mp_t> mp; // array of material properties
      extern std::vector<double> stt_rj; // array of adiabatic spin torques
      extern std::vector<double> stt_pj; // array of non-adiabatic spin torques
      extern std::vector<double> stt_polarization_unit_vector; // stt spin polarization direction
      extern std::vector<double> sot_rj; // array of adiabatic spin torques
      extern std::vector<double> sot_pj; // array of non-adiabatic spin torques
      extern std::vector<double> sot_polarization_unit_vector; // sot spin polarization direction

   } // end of internal namespace
} // end of sim namespace

#endif //SIM_INTERNAL_H_
