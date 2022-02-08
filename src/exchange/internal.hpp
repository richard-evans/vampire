//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef EXCHANGE_INTERNAL_H_
#define EXCHANGE_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the exchange module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "exchange.hpp"

// exchange module headers
#include "internal.hpp"

namespace exchange{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------
      //-----------------------------------------------------------------------------
      // materials class for storing exchange material parameters
      //-----------------------------------------------------------------------------
      class mp_t{

         private:

         public:

            // variables
            std::vector<double> bqe; // Biquadratic exchange constant
            std::vector<double> dmi; // Dzyaloshinskii-Moriya interaction constant

            // constructor
            mp_t (const unsigned int max_materials = 100)
            {
               // resize arrays to correct size
               bqe.resize(max_materials, 0.0); // initialise biquadratic constants to zero
               dmi.resize(max_materials, 0.0); // initialise pair anisotropy constants to zero

            }; // end of constructor

      }; // end of exchange::internal::mp class

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------
      extern std::vector<internal::mp_t> mp; // array of material properties

      extern bool enable_dmi; // flag to enable dmi calculation

      extern double dmi_cutoff_range; // cutoff range for DMI calculation (Ångstroms)

      extern exchange_t exchange_type; // exchange type to use in simulation
      extern exchange_t biquadratic_exchange_type; // biquadratic exchange type to use in simulation

      extern exchange_t minimum_needed_exchange_type; // minimum support required for bilinear exchange type (for vectorial constants and DMI)

      extern bool use_material_exchange_constants; // flag to enable material exchange parameters
      extern bool use_material_biquadratic_exchange_constants; // flag to enable material exchange parameters

      // simple vector and tensor class definitions
      class value_t{
      	public:
      	double Jij;

      	// constructor
      	value_t():
      		Jij(0.0)
      	{
      	};
      };

      class vector_t{
      	public:
      	double Jij[3];

      	// constructor
      	vector_t()
      	{
      		Jij[0] = 0.0;
      		Jij[1] = 0.0;
      		Jij[2] = 0.0;
      	};
      };

      class tensor_t{
      	public:
      	double Jij[3][3];

      	// constructor
      	tensor_t()
      	{
      		Jij[0][0] = 0.0;
      		Jij[0][1] = 0.0;
      		Jij[0][2] = 0.0;

      		Jij[1][0] = 0.0;
      		Jij[1][1] = 0.0;
      		Jij[1][2] = 0.0;

      		Jij[2][0] = 0.0;
      		Jij[2][1] = 0.0;
      		Jij[2][2] = 0.0;
      	};
      };

      extern std::vector <int> biquadratic_neighbour_list_array; // 1D list of biquadratic neighbours
      extern std::vector <int> biquadratic_neighbour_interaction_type_array; // 1D list of biquadratic exchange interaction types
      extern std::vector <int> biquadratic_neighbour_list_start_index; // list of first biquadratic neighbour for atom i
      extern std::vector <int> biquadratic_neighbour_list_end_index;   // list of last biquadratic neighbour for atom i

      extern std::vector <exchange::internal::value_t > bq_i_exchange_list; // list of isotropic biquadratic exchange constants
      extern std::vector <exchange::internal::vector_t> bq_v_exchange_list; // list of vectorial biquadratic exchange constants
      extern std::vector <exchange::internal::tensor_t> bq_t_exchange_list; // list of tensorial biquadratic exchange constants

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------
      void calculate_dmi(std::vector<std::vector <neighbours::neighbour_t> >& cneighbourlist);
      void unroll_exchange_interactions();
      void unroll_normalised_exchange_interactions();
      void unroll_normalised_biquadratic_exchange_interactions();
      void exchange_fields(const int start_index, // first atom for exchange interactions to be calculated
                           const int end_index, // last +1 atom to be calculated
                           const std::vector<int>& neighbour_list_start_index,
                           const std::vector<int>& neighbour_list_end_index,
                           const std::vector<int>& type_array, // type for atom
                           const std::vector<int>& neighbour_list_array, // list of interactions between atoms
                           const std::vector<int>& neighbour_interaction_type_array, // list of interaction type for each pair of atoms with value given in exchange list
                           const std::vector<zval_t>& i_exchange_list, // list of isotropic exchange constants
                           const std::vector<zvec_t>& v_exchange_list, // list of vectorial exchange constants
                           const std::vector<zten_t>& t_exchange_list, // list of tensorial exchange constants
                           const std::vector<double>& spin_array_x, // spin vectors for atoms
                           const std::vector<double>& spin_array_y,
                           const std::vector<double>& spin_array_z,
                           std::vector<double>& field_array_x, // field vectors for atoms
                           std::vector<double>& field_array_y,
                           std::vector<double>& field_array_z);

      void biquadratic_exchange_fields(const int start_index, // first atom for exchange interactions to be calculated
                                       const int end_index, // last +1 atom to be calculated
                                       const std::vector<int>& neighbour_list_start_index,
                                       const std::vector<int>& neighbour_list_end_index,
                                       const std::vector<int>& type_array, // type for atom
                                       const std::vector<int>& neighbour_list_array, // list of interactions between atoms
                                       const std::vector<int>& neighbour_interaction_type_array, // list of interaction type for each pair of atoms with value given in exchange list
                                       const std::vector<value_t>&  bq_i_exchange_list, // list of isotropic biquadratic exchange constants
                                       const std::vector<vector_t>& bq_v_exchange_list, // list of vectorial biquadratic exchange constants
                                       const std::vector<tensor_t>& bq_t_exchange_list, // list of tensorial biquadratic exchange constants
                                       const std::vector<double>& spin_array_x, // spin vectors for atoms
                                       const std::vector<double>& spin_array_y,
                                       const std::vector<double>& spin_array_z,
                                       std::vector<double>& field_array_x, // field vectors for atoms
                                       std::vector<double>& field_array_y,
                                       std::vector<double>& field_array_z);

      void initialize_biquadratic_exchange();

   } // end of internal namespace

} // end of exchange namespace

#endif //EXCHANGE_INTERNAL_H_
