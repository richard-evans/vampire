//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Fried-Conrad Weber 2025. All rights reserved.
//
//   Email: fried-conrad.weber@uni-potsdam.de
//
//------------------------------------------------------------------------------
//

#ifndef QUANTUM_INTERNAL_H_
#define QUANTUM_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the quantum module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <vector>

// Vampire headers
#include "quantum.hpp"

// quantum module headers
#include "internal.hpp"

namespace quantum{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      //-----------------------------------------------------------------------------
      // internal materials class for storing material parameters
      //-----------------------------------------------------------------------------
      class mp_t{

          private:
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
                   }
                   // getting function
                   double get(){
                      return value;
                   }
                   // checking function
                   bool is_set(){
                      return setf;
                   }
             };

          public:

             //------------------------------
             // material parameter variables
             //------------------------------
             set_double_t A;
             set_double_t Gamma;
             set_double_t omega0;
             set_double_t S0;
      };

      //------------------------------------------------------------------------
      // Shared variables inside quantum module
      //------------------------------------------------------------------------

      extern bool enabled; // bool to enable module

      extern std::vector<internal::mp_t> mp; // array of material properties

      extern noise_t noise_type; // Default to Quantum noise
      extern int window_size; // Default window size (initially set to zero, later initialised to simulation length)
      extern int M_decimation; // No decimation (interpolation in German -> English translation)

      extern std::vector<double> material_A_array; // compact array of Lorentzian amplitude for each material
      extern std::vector<double> material_gamma_array; // compact array of Lorentzian gamma for each material
      extern std::vector<double> material_omega0_array; // compact array of Lorentzian omega0 for each material

      extern std::vector<double> coarse_noise_field; // final generated noise (1D list of noise for all atoms and dimensions)

      extern std::vector<std::size_t> atom_idx_x; // indexing for accessing noise for each atom and spatial dimension in coarse array
      extern std::vector<std::size_t> atom_idx_y;
      extern std::vector<std::size_t> atom_idx_z;

      // Integration arrays for auxiliary variables
      extern std::vector<double> x_w_array;
      extern std::vector<double> y_w_array;
      extern std::vector<double> z_w_array;
      extern std::vector<double> x_v_array;
      extern std::vector<double> y_v_array;
      extern std::vector<double> z_v_array;

      // Integration arrays for RK4 integration
      extern std::vector<std::vector<double>> k1_storage;
      extern std::vector<std::vector<double>> k2_storage;
      extern std::vector<std::vector<double>> k3_storage;
      extern std::vector<std::vector<double>> k4_storage;
      extern std::vector<std::vector<double>> y_pred_storage;
      extern std::vector<std::vector<double>> y_in_storage;

      //------------------------------------------------------------------------
      // Function prototypes
      //------------------------------------------------------------------------
      double PSD(const double& omega, const double& T);
      double estimate_cutoff_omega_cdf(double T, double target_frac);
      void calculate_random_fields_windowed(int realizations, int n_fine, double dt_fine, int M, double T, int n_coarse_total, std::vector<double>& noise_field);
      double get_noise(const std::vector<double>& coarse_noise, double fine_step_idx, int M, size_t atom_idx);
      void assign_unique_indices(int n_coarse);
             double test;

             // constructor
             mp_t (const unsigned int max_materials = 100):
                test(0.0) // constructor initialisation of test variable
             {
                // constructor body for initialising more complex data/arrays
             }; // end of constructor

       }; // end of internal::mp class

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------

      extern bool enabled; // bool to enable module

      extern std::vector<internal::mp_t> mp; // array of material properties

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

   } // end of internal namespace

} // end of quantum namespace

#endif //QUANTUM_INTERNAL_H_
