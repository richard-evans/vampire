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
         set_double_t slonczewski_aj; // adiabatic torque
         set_double_t slonczewski_bj; // field torque (non-adiabatic)
         set_double_t SOT_DL; // damping-like SOT torque
         set_double_t SOT_FL; // field-like SOT torque (non-adiabatic)
      };

      extern std::vector<sim::internal::mp_t> mp; // array of material properties
      extern std::vector<double> slonczewski_aj; // array of adiabatic spin torques
      extern std::vector<double> slonczewski_bj; // array of non-adiabatic spin torques
      extern std::vector<double> slonczewski_spin_polarization_unit_vector; // spin polarization direction
      extern std::vector<double> SOT_DL; // array of adiabatic spin torques
      extern std::vector<double> SOT_FL; // array of non-adiabatic spin torques
      extern std::vector<double> SOT_spin_polarization_unit_vector; // spin polarization direction


   } // end of internal namespace
} // end of sim namespace

#endif //SIM_INTERNAL_H_
