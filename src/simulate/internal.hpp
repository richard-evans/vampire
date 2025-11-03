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
      // Internal data types used for simulation module
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
         set_double_t stt_asm; // spin tranfer torque asymmetry
         set_double_t stt_rj;  // spin tranfer relaxation torque
         set_double_t stt_pj;  // spin transfer precession torque
         set_double_t sot_asm; // spin orbit torque asymmetry
         set_double_t sot_rj;  // spin orbit relaxation torque
         set_double_t sot_pj;  // spin orbit precession torque
         set_double_t vcmak;   // voltage controlled anisotropy coefficient
         set_double_t lsf_second_order_coefficient; // Second order LSF coefficient
         set_double_t lsf_fourth_order_coefficient; // Fourth order LSF coefficient
         set_double_t lsf_sixth_order_coefficient; // Sixth order LSF coefficient
         set_double_t A; // Quantum thermostat parameter
         set_double_t Gamma; // Quantum thermostat parameter
         set_double_t omega0; // Quantum thermostat parameter
         set_double_t S0; // Quantum thermostat parameter
      };

      //-----------------------------------------------------------------------------
      // Internal shared variables used for the simulation
      //-----------------------------------------------------------------------------
      extern bool enable_spin_torque_fields; // flag to enable spin torque fields
      extern bool enable_vcma_fields;        // flag to enable voltage-controlled anisotropy fields

      extern std::vector<sim::internal::mp_t> mp; // array of material properties

      extern std::vector<double> stt_asm; // array of spin transfer torque asymmetry
      extern std::vector<double> stt_rj; // array of adiabatic spin torques
      extern std::vector<double> stt_pj; // array of non-adiabatic spin torques
      extern std::vector<double> stt_polarization_unit_vector; // stt spin polarization direction

      extern std::vector<double> sot_asm; // array of spin orbit torque asymmetry
      extern std::vector<double> sot_rj; // array of adiabatic spin torques
      extern std::vector<double> sot_pj; // array of non-adiabatic spin torques
      extern std::vector<double> sot_polarization_unit_vector; // sot spin polarization direction

      extern std::vector<double> vcmak;   // voltage controlled anisotropy coefficient

      extern std::vector<double> lsf_second_order_coefficient;
      extern std::vector<double> lsf_fourth_order_coefficient; // LSF coefficients
      extern std::vector<double> lsf_sixth_order_coefficient;

      // shared Functions
      void llg_quantum_step();

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------
      extern void initialize_modules();
      extern void increment_time();
      extern void lsf_step();
      extern void lsf_rk4_step();

   } // end of internal namespace

   namespace LLGQ_arrays{

   //---------------------------------------------------------------------------
   // Local arrays and variables
   //---------------------------------------------------------------------------

   // Auxiallary variables
   extern std::vector <double> x_w_array;
   extern std::vector <double> y_w_array;
   extern std::vector <double> z_w_array;
   extern std::vector <double> x_v_array;
   extern std::vector <double> y_v_array;
   extern std::vector <double> z_v_array;

   // Storage arrays for RK4
   extern std::vector<std::vector<double>> k1_storage;
   extern std::vector<std::vector<double>> k2_storage;
   extern std::vector<std::vector<double>> k3_storage;
   extern std::vector<std::vector<double>> k4_storage;
   extern std::vector<std::vector<double>> y_pred_storage;
   extern std::vector<std::vector<double>> y_in_storage;

   // Arrays for noise generation (now storing coarse-grained noise)
   extern std::vector<double> coarse_noise_field;
   extern std::vector<double> sqrt_PSD_buffer; // Restored buffer
   extern double noise_index;
   extern int M_decimation;    //< Decimation factor for noise interpolation

   // Indices for random fields
   extern std::vector<double> atom_idx_x;
   extern std::vector<double> atom_idx_y;
   extern std::vector<double> atom_idx_z;

   extern bool LLG_set; ///< Flag to define state of LLG arrays (initialised/uninitialised)

   }
} // end of sim namespace

#endif //SIM_INTERNAL_H_
