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

namespace quantum{

   namespace internal{
      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------
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
      //-----------------------------------------------------------------------------
      // materials class for storing exchange material parameters
      //-----------------------------------------------------------------------------
      class mp_t{

         private:

         public:

            //----------------
            // variables
            //----------------
            set_double_t gamma;      // Lorentzian width parameter
            set_double_t omega0;     // Lorentzian central frequency parameter

            // constructor
            mp_t (const unsigned int max_materials = 100)
            {
               gamma.set(0.0); // default value
               omega0.set(0.0); // default value
            }; // end of constructor

      };

      //------------------------------------------------------------------------
      // Shared variables inside quantum module
      //------------------------------------------------------------------------

      enum noise_t {
         classical,
         quantum_zero,
         quantum_no_zero
      };

      extern bool initialised; // bool to enable module

      extern std::vector<internal::mp_t> mp; // array of material properties

      extern noise_t noise_type; // Default to Quantum noise
      extern uint64_t window_size; // Default window size (initially set to zero, later initialised to simulation length)
      extern uint64_t M_decimation; // No decimation (interpolation in German -> English translation)

      extern std::vector<double> material_A_array; // compact array of Lorentzian amplitude for each material
      extern std::vector<double> material_gamma_array; // compact array of Lorentzian gamma for each material
      extern std::vector<double> material_omega0_array; // compact array of Lorentzian omega0 for each material

      extern std::vector<double> coarse_noise_field; // final generated noise (1D list of noise for all atoms and dimensions)

      extern std::vector<std::size_t> atom_idx_x; // indexing for accessing noise for each atom and spatial dimension in coarse array
      extern std::vector<std::size_t> atom_idx_y;
      extern std::vector<std::size_t> atom_idx_z;

      extern int noise_index; // index for accessing noise arrays (at each fine time step)

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

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

      double PSD(const double& omega,
                 const double& T);

      void calculate_noise(int realizations,
                           int n_fine,
                           double dt_fine,
                           int M,
                           double T,
                           int n_coarse_total,
                           std::vector<double>& noise_field);


      double get_noise(const std::vector<double>& coarse_noise,
                       double fine_step_idx,
                       int M,
                       size_t atom_idx);

      void assign_unique_indices(int n_coarse);

   } // end of internal namespace

} // end of quantum namespace

#endif //QUANTUM_INTERNAL_H_
