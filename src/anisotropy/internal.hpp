//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sam Westmoreland and Richard Evans 2017. All rights reserved.
//
//   Email: sw766@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef ANISOTROPY_INTERNAL_H_
#define ANISOTROPY_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the anisotropy module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <vector>

// Vampire headers
#include "material.hpp"
#include "anisotropy.hpp"

// anisotropy module headers
#include "internal.hpp"

namespace anisotropy{

   namespace internal{

      //-----------------------------------------------------------------------------
      // Internal data type definitions
      //-----------------------------------------------------------------------------

      //--------------------------------------------------------------------
      // Class to contain parameters for lattice anisotropy calculation.
      //
      // Tabulated values are read from a file and added point-wise to
      // the class. During variable initialisation interpolating functions
      // are determined to calculate k(T)
      //
      class lattice_anis_t{

         private:

            unsigned int Tmax; // maximum array value in tabulated function
            double k_Tmax; // value of anisotropy at k_Tmax (used for all T > Tmax)

            std::vector<unsigned int> T; // input temperature points
            std::vector<double> k; // input lattice anisotropy points
            std::vector<double> m; // calculated m value
            std::vector<double> c; // calculated c value

         public:

            void add_point(double, double);
            void set_interpolation_table();
            double get_lattice_anisotropy_constant(double);
            void output_interpolated_function(int);

      };

      //-----------------------------------------------------------------------------
      // materials class for storing anisotropy material parameters
      //-----------------------------------------------------------------------------
      class mp_t{

         private:

         public:

            // variables
            double ku2; // second order uniaxial anisotropy constant (Ku1)
            double ku4; // fourth order uniaxial anisotropy constant (Ku2)
            double ku6; // sixth order uniaxial anisotropy constant  (Ku3)

            double kc4; // fourth order cubic anisotropy constant (Kc1)
            double kc6; // sixth order cubic anisotropy constant (Kc2)

            double k_lattice; // uniaxial lattice anisotropy constant

            std::vector<double> kij; // surface/Neel anisotropy pair constant

            std::vector<double> ku_vector; // unit vector defining axis for uniaxial anisotropy
            std::vector<double> kc_vector; // unit vector defining axis for cubic anisotropy

            std::vector<double> ku_tensor; // uniaxial second order anisotropy tensor
            std::vector<double> kc_tensor; // cubic fourth order anisotropy tensor

            lattice_anis_t lattice_anisotropy; // class containing lattice anisotropy data

            bool random_anisotropy; // flag to control random anisotropy by material
      		bool random_grain_anisotropy; // flag to control random anisotropy by grain

            // constructor
            mp_t (const unsigned int max_materials = 100):
            	ku2(0.0), // set initial value of ku2 to zero
               ku4(0.0), // set initial value of ku4 to zero
               ku6(0.0), // set initial value of ku6 to zero
               kc4(0.0), // set initial value of kc4 to zero
               kc6(0.0), // set initial value of kc6 to zero
               random_anisotropy(false), // disable random anisotropy
               random_grain_anisotropy(false) // disable random grain anisotropy
            {
               // resize arrays to correct size
               kij.resize(max_materials, 0.0); // initialise pair anisotropy constants to zero

               // set default uniaxial and cubic directions
               ku_vector.resize(3); // resize to three elements

               ku_vector[0] = 0.0; // set direction alon [0,0,1]
               ku_vector[1] = 0.0;
               ku_vector[2] = 1.0;

               // set default uniaxial and cubic directions
               kc_vector.resize(3); // resize to three elements

               kc_vector[0] = 0.0; // set direction alon [0,0,1]
               kc_vector[1] = 0.0;
               kc_vector[2] = 1.0;

               // set tensors as empty by default
               ku_tensor.resize(9, 0.0);
               kc_tensor.resize(9, 0.0);

            }; // end of constructor

      }; // end of anisotropy::internal::mp class

      //-----------------------------------------------------------------------------
      // Internal shared variables used for creation
      //-----------------------------------------------------------------------------

      extern std::vector<internal::mp_t> mp; // array of material properties

      extern bool initialised; // check module has been initialised

      extern bool enable_second_order_tensor; // Flag to enable calculation of second order tensor anisotropy
      extern bool enable_fourth_order_tensor; // Flag to enable calculation of second order tensor anisotropy
      extern bool enable_sixth_order_tensor; // Flag to enable calculation of second order tensor anisotropy

      extern bool enable_neel_anisotropy; // Flag to turn on Neel anisotropy calculation (memory intensive at startup)
      extern bool enable_lattice_anisotropy; // Flag to turn on lattice anisotropy calculation
      extern bool enable_random_anisotropy; // Flag to enable random anisitropy initialisation

      // arrays for storing 1D collapsed tensors
      extern std::vector<double> second_order_tensor;
      extern std::vector<double> fourth_order_tensor;
      extern std::vector<double> sixth_order_tensor;

      extern bool native_neel_anisotropy_threshold; // enables site-dependent surface threshold
      extern unsigned int neel_anisotropy_threshold; // global threshold for surface atoms
      extern double nearest_neighbour_distance; // Control surface anisotropy nearest neighbour distance

      //-------------------------------------------------------------------------
      // internal function declarations
      //-------------------------------------------------------------------------
      extern void uniaxial_second_order(const unsigned int num_atoms, std::vector<int>& atom_material_array, std::vector<double>& inverse_mu_s);
      extern void uniaxial_fourth_order(const unsigned int num_atoms, std::vector<int>& atom_material_array, std::vector<double>& inverse_mu_s);

   } // end of internal namespace

   //-------------------------------------------------------------------------
   // function declarations
   //-------------------------------------------------------------------------

   //extern int calculate_energies();
   //extern double second_order_tensor_anisotropy();
   //extern double third_order_tensor_anisotropy();
   //extern double fourth_order_tensor_anisotropy();

   //-------------------------------------------------------------------------
   // simple inline function to convert atom,i,j into 1D tensor coordinates
   //-------------------------------------------------------------------------
   inline unsigned int index(const unsigned int atom, const unsigned int i, const unsigned int j){
      return 9*atom + 3*i + j;
   }

} // end of anisotropy namespace

#endif //ANISOTROPY_INTERNAL_H_
