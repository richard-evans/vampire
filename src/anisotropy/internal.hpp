//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Roberto Moreno Ortega, Sam Westmoreland and Richard Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
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

      //---------------------------------
      // struct for storing unit vectors
      //---------------------------------
      struct evec_t{
         double x;
         double y;
         double z;
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

            double k4r; // fourth order rotational anisotropy constant (k4r)

            double k_lattice; // uniaxial lattice anisotropy constant

            std::vector<double> kij; // surface/Neel anisotropy pair constant

            std::vector<double> ku_vector; // unit vector defining axis for uniaxial anisotropy

            std::vector<double> u1_vector; // unit vector defining axis for uniaxial anisotropy
            std::vector<double> u2_vector; // unit vector defining axis for uniaxial anisotropy

            std::vector<double> kc_vector1; // first unit vector defining axis for cubic anisotropy
            std::vector<double> kc_vector2; // second unit vector defining axis for cubic anisotropy
            std::vector<double> kc_vector3; // third unit vector defining axis for cubic anisotropy

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
               k4r(0.0), // set initial value of k4r to zero
               k_lattice(0.0), // set initial value of k_lattice to zero
               random_anisotropy(false), // disable random anisotropy
               random_grain_anisotropy(false) // disable random grain anisotropy
            {
               // resize arrays to correct size
               kij.resize(max_materials, 0.0); // initialise pair anisotropy constants to zero

               // set default uniaxial and cubic directions
               ku_vector.resize(3); // resize to three elements

               ku_vector[0] = 0.0; // set direction along [0,0,1]
               ku_vector[1] = 0.0;
               ku_vector[2] = 1.0;

               const double oneosqrt2 = 1.0/sqrt(2.0);

               u1_vector.resize(3); // resize to three elements
               u2_vector.resize(3); // resize to three elements

               u1_vector[0] = oneosqrt2*1.0; // set direction along [1,1,0]
               u1_vector[1] = oneosqrt2*1.0;
               u1_vector[2] = 0.0;

               u2_vector[0] = oneosqrt2*1.0; // set direction along [1,-1,0]
               u2_vector[1] = -oneosqrt2*1.0;
               u2_vector[2] = 0.0;

               // set default uniaxial and cubic directions
               kc_vector1.resize(3); // resize to three elements
               kc_vector2.resize(3); // resize to three elements
               kc_vector3.resize(3); // resize to three elements

               kc_vector1[0] = 1.0; // set direction alon [1,0,0]
               kc_vector1[1] = 0.0;
               kc_vector1[2] = 0.0;

               kc_vector2[0] = 0.0; // set direction alon [0,1,0]
               kc_vector2[1] = 1.0;
               kc_vector2[2] = 0.0;

               kc_vector3[0] = 0.0; // set direction alon [0,0,1]
               kc_vector3[1] = 0.0;
               kc_vector3[2] = 1.0;

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

      extern bool enable_uniaxial_second_order; // Flag to enable calculation of second order uniaxial anisotropy
      extern bool enable_uniaxial_fourth_order; // Flag to enable calculation of fourth order uniaxial anisotropy
      extern bool enable_biaxial_fourth_order_simple; // Flag to enable calculation of simplified fourth order biaxial anisotropy
      extern bool enable_uniaxial_sixth_order;  // Flag to enable calculation of sixth order uniaxial anisotropy

      extern bool enable_cubic_fourth_order;    // Flag to enable calculation of fourth order cubic anisotropy
      extern bool enable_cubic_sixth_order;     // Flag to enable calculation of sixth order cubic  anisotropy
      extern bool enable_cubic_fourth_order_rotation; // Flag to enable calculation of rotated cubic anisotropy

      extern bool enable_fourth_order_rotational; // Flag to enable 4th order rotational anisotropy

      extern bool enable_triaxial_anisotropy_rotated;
      extern bool enable_triaxial_fourth_order_rotated;
      extern bool enable_triaxial_anisotropy;
      extern bool enable_triaxial_fourth_order;

      extern bool enable_neel_anisotropy; // Flag to turn on Neel anisotropy calculation (memory intensive at startup)
      extern bool enable_lattice_anisotropy; // Flag to turn on lattice anisotropy calculation
      extern bool enable_random_anisotropy; // Flag to enable random anisitropy initialisation

      // arrays for storing 1D collapsed Neel tensor
      extern std::vector<double> neel_tensor;

      // arrays for storing unrolled anisotropy constants in Tesla
      extern std::vector<double> ku2;
      extern std::vector<double> ku4;
      extern std::vector<double> ku6;
      extern std::vector<double> kc4;
      extern std::vector<double> kc6;
      extern std::vector<double> k4r;

      // unrolled arrays for storing easy axes for each material
      extern std::vector<evec_t> ku_vector; // 001 easy axis direction

      extern bool native_neel_anisotropy_threshold;  // enables site-dependent surface threshold
      extern unsigned int neel_anisotropy_threshold; // global threshold for surface atoms
      extern double nearest_neighbour_distance;      // Control surface anisotropy nearest neighbour distance
      extern bool neel_range_dependent;              // Enable range dependent Neel anisotropy Lij = L0 exp(-F(r-r0)/r0)
      extern double neel_exponential_range;          // r0 value for range dependence of Neel anisotropy
      extern double neel_exponential_factor;         // F value for range dependence of Neel anisotropy

      extern std::vector<bool> triaxial_second_order_fixed_basis;
      extern std::vector<bool> triaxial_fourth_order_fixed_basis;

      extern std::vector<double> ku_triaxial_vector_x; // unit vector defining axis for uniaxial anisotropy
      extern std::vector<double> ku_triaxial_vector_y; // unit vector defining axis for uniaxial anisotropy
      extern std::vector<double> ku_triaxial_vector_z; // unit vector defining axis for uniaxial anisotropy

      extern std::vector<double> ku4_triaxial_vector_x; // unit vector defining axis for uniaxial anisotropy
      extern std::vector<double> ku4_triaxial_vector_y; // unit vector defining axis for uniaxial anisotropy
      extern std::vector<double> ku4_triaxial_vector_z; // unit vector defining axis for uniaxial anisotropy

      //basis vectors for second order triaxial - must be orthogonality
      extern std::vector < double > ku_triaxial_basis1x;
      extern std::vector < double > ku_triaxial_basis1y;
      extern std::vector < double > ku_triaxial_basis1z;

      extern std::vector < double > ku_triaxial_basis2x;
      extern std::vector < double > ku_triaxial_basis2y;
      extern std::vector < double > ku_triaxial_basis2z;

      extern std::vector < double > ku_triaxial_basis3x;
      extern std::vector < double > ku_triaxial_basis3y;
      extern std::vector < double > ku_triaxial_basis3z;

      //basis vectors for fourth order triaxial - must be orthogonality
      extern std::vector < double > ku4_triaxial_basis1x;
      extern std::vector < double > ku4_triaxial_basis1y;
      extern std::vector < double > ku4_triaxial_basis1z;

      extern std::vector < double > ku4_triaxial_basis2x;
      extern std::vector < double > ku4_triaxial_basis2y;
      extern std::vector < double > ku4_triaxial_basis2z;

      extern std::vector < double > ku4_triaxial_basis3x;
      extern std::vector < double > ku4_triaxial_basis3y;
      extern std::vector < double > ku4_triaxial_basis3z;

      // arrays for storing unrolled parameters for lattice anisotropy
      extern std::vector< double > klattice; // anisotropy constant
      extern std::vector< double > klattice_array; // array for unrolled anisotropy including temperature dependence

      //-------------------------------------------------------------------------
      // internal function declarations
      //-------------------------------------------------------------------------
      void uniaxial_second_order_fields(std::vector<double>& spin_array_x,
                                        std::vector<double>& spin_array_y,
                                        std::vector<double>& spin_array_z,
                                        std::vector<int>&    atom_material_array,
                                        std::vector<double>& field_array_x,
                                        std::vector<double>& field_array_y,
                                        std::vector<double>& field_array_z,
                                        const int start_index,
                                        const int end_index);


      void triaxial_second_order_fields_fixed_basis(std::vector<double>& spin_array_x,
                                                     std::vector<double>& spin_array_y,
                                                     std::vector<double>& spin_array_z,
                                                     std::vector<int>&    atom_material_array,
                                                     std::vector<double>& field_array_x,
                                                     std::vector<double>& field_array_y,
                                                     std::vector<double>& field_array_z,
                                                     const int start_index,
                                                     const int end_index);

      void triaxial_second_order_fields(std::vector<double>& spin_array_x,
                                                    std::vector<double>& spin_array_y,
                                                    std::vector<double>& spin_array_z,
                                                    std::vector<int>&    atom_material_array,
                                                    std::vector<double>& field_array_x,
                                                    std::vector<double>& field_array_y,
                                                    std::vector<double>& field_array_z,
                                                    const int start_index,
                                                    const int end_index);

      void triaxial_fourth_order_fields_fixed_basis(std::vector<double>& spin_array_x,
                                                  std::vector<double>& spin_array_y,
                                                  std::vector<double>& spin_array_z,
                                                  std::vector<int>&    atom_material_array,
                                                  std::vector<double>& field_array_x,
                                                  std::vector<double>& field_array_y,
                                                  std::vector<double>& field_array_z,
                                                  const int start_index,
                                                  const int end_index);

      void triaxial_fourth_order_fields(std::vector<double>& spin_array_x,
                                        std::vector<double>& spin_array_y,
                                        std::vector<double>& spin_array_z,
                                        std::vector<int>&    atom_material_array,
                                        std::vector<double>& field_array_x,
                                        std::vector<double>& field_array_y,
                                        std::vector<double>& field_array_z,
                                        const int start_index,
                                        const int end_index);

      void uniaxial_fourth_order_fields(std::vector<double>& spin_array_x,
                                        std::vector<double>& spin_array_y,
                                        std::vector<double>& spin_array_z,
                                        std::vector<int>&    atom_material_array,
                                        std::vector<double>& field_array_x,
                                        std::vector<double>& field_array_y,
                                        std::vector<double>& field_array_z,
                                        const int start_index,
                                        const int end_index);

      void biaxial_fourth_order_simple_fields(std::vector<double>& spin_array_x,
                                        std::vector<double>& spin_array_y,
                                        std::vector<double>& spin_array_z,
                                        std::vector<int>&    atom_material_array,
                                        std::vector<double>& field_array_x,
                                        std::vector<double>& field_array_y,
                                        std::vector<double>& field_array_z,
                                        const int start_index,
                                        const int end_index);

      void uniaxial_sixth_order_fields( std::vector<double>& spin_array_x,
                                        std::vector<double>& spin_array_y,
                                        std::vector<double>& spin_array_z,
                                        std::vector<int>&    atom_material_array,
                                        std::vector<double>& field_array_x,
                                        std::vector<double>& field_array_y,
                                        std::vector<double>& field_array_z,
                                        const int start_index,
                                        const int end_index);

      void rotational_fourth_order_fields_fixed_basis( std::vector<double>& spin_array_x,
                                        std::vector<double>& spin_array_y,
                                        std::vector<double>& spin_array_z,
                                        std::vector<int>&    atom_material_array,
                                        std::vector<double>& field_array_x,
                                        std::vector<double>& field_array_y,
                                        std::vector<double>& field_array_z,
                                        const int start_index,
                                        const int end_index);

      void cubic_fourth_order_fields(std::vector<double>& spin_array_x,
                                     std::vector<double>& spin_array_y,
                                     std::vector<double>& spin_array_z,
                                     std::vector<int>&    atom_material_array,
                                     std::vector<double>& field_array_x,
                                     std::vector<double>& field_array_y,
                                     std::vector<double>& field_array_z,
                                     const int start_index,
                                     const int end_index);

      void cubic_fourth_order_rotation_fields(std::vector<double>& spin_array_x,
                                              std::vector<double>& spin_array_y,
                                              std::vector<double>& spin_array_z,
                                              std::vector<int>&    atom_material_array,
                                              std::vector<double>& field_array_x,
                                              std::vector<double>& field_array_y,
                                              std::vector<double>& field_array_z,
                                              const int start_index,
                                              const int end_index);

      void cubic_sixth_order_fields( std::vector<double>& spin_array_x,
                                     std::vector<double>& spin_array_y,
                                     std::vector<double>& spin_array_z,
                                     std::vector<int>&    atom_material_array,
                                     std::vector<double>& field_array_x,
                                     std::vector<double>& field_array_y,
                                     std::vector<double>& field_array_z,
                                     const int start_index,
                                     const int end_index);

      void neel_fields( std::vector<double>& spin_array_x,
                        std::vector<double>& spin_array_y,
                        std::vector<double>& spin_array_z,
                        std::vector<int>&    atom_material_array,
                        std::vector<double>& field_array_x,
                        std::vector<double>& field_array_y,
                        std::vector<double>& field_array_z,
                        const int start_index,
                        const int end_index);

      void lattice_fields(std::vector<double>& spin_array_x,
                          std::vector<double>& spin_array_y,
                          std::vector<double>& spin_array_z,
                          std::vector<int>&    type_array,
                          std::vector<double>& field_array_x,
                          std::vector<double>& field_array_y,
                          std::vector<double>& field_array_z,
                          const int start_index,
                          const int end_index,
                          const double temperature);

      double uniaxial_second_order_energy( const int atom,
                                           const int mat,
                                           const double sx,
                                           const double sy,
                                           const double sz);

      double triaxial_second_order_energy_fixed_basis(const int atom,
                                                      const int mat,
                                                      const double sx,
                                                      const double sy,
                                                      const double sz);

      double triaxial_fourth_order_energy_fixed_basis(const int atom,
                                                      const int mat,
                                                      const double sx,
                                                      const double sy,
                                                      const double sz);

      double triaxial_second_order_energy(const int atom,
                                          const int mat,
                                          const double sx,
                                          const double sy,
                                          const double sz);

      double triaxial_fourth_order_energy(const int atom,
                                          const int mat,
                                          const double sx,
                                          const double sy,
                                          const double sz);


      double uniaxial_fourth_order_energy(const int atom,
                                          const int mat,
                                          const double sx,
                                          const double sy,
                                          const double sz);

      double biaxial_fourth_order_simple_energy(const int atom,
                                          const int mat,
                                          const double sx,
                                          const double sy,
                                          const double sz);

      double uniaxial_sixth_order_energy( const int atom,
                                          const int mat,
                                          const double sx,
                                          const double sy,
                                          const double sz);

      double rotational_fourth_order_energy_fixed_basis( const int atom,
                                                         const int mat,
                                                         const double sx,
                                                         const double sy,
                                                         const double sz);

      double cubic_fourth_order_energy(const int atom,
                                       const int mat,
                                       const double sx,
                                       const double sy,
                                       const double sz);

      double cubic_fourth_order_rotation_energy(const int atom,
                                       const int mat,
                                       const double sx,
                                       const double sy,
                                       const double sz);


      double cubic_sixth_order_energy( const int atom,
                                       const int mat,
                                       const double sx,
                                       const double sy,
                                       const double sz);

      double neel_energy( const int atom,
                          const int mat,
                          const double sx,
                          const double sy,
                          const double sz);

      double lattice_energy(const int atom, const int mat, const double sx, const double sy, const double sz, const double temperature);

      void initialise_neel_anisotropy_tensor(std::vector <std::vector <bool> >& nearest_neighbour_interactions_list,
                                             std::vector<std::vector <neighbours::neighbour_t> >& cneighbourlist);

   } // end of internal namespace

   //-------------------------------------------------------------------------
   // function declarations
   //-------------------------------------------------------------------------

   //-------------------------------------------------------------------------
   // simple inline function to convert atom,i,j into 1D tensor coordinates
   //-------------------------------------------------------------------------
   inline unsigned int index(const unsigned int atom, const unsigned int i, const unsigned int j){
      return 9*atom + 3*i + j;
   }

} // end of anisotropy namespace

#endif //ANISOTROPY_INTERNAL_H_
