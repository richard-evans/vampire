//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2020. All rights reserved.
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
#include "errors.hpp"
#include "vio.hpp"

// exchange module headers
#include "internal.hpp"

namespace exchange{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------
      // data type for storing long range user definable long range exchange
      // interactions. In general data can be num_materials * num_materials * 10 * 9
      // elements large, and so we store the data in compact arrays so that only
      // defined materials are stored, and only the righ number of values.
      // access is horrible so class includes set and get functions to make things
      // easier programatically. In the end these values are just unrolled into the
      // usual exchange list.
      //
      class exchange_matrix_4D_t{

         private:
            int max_materials; // number of materials needed to store user input data
            int max_neighbours; // number of neighbour needed (determines range of interaction)
            // lovely 4D array for storing [materiali][materialj][neighbour][Jijs]
            std::vector< std::vector< std::vector< std::vector <double> > > > exchange_matrix_array;

            // function to resize storage array
            void resize(const int num_materials, const int num_neighbours){

               //std::cout << "Setting storage size to " << num_materials << " materials and " << num_neighbours << " neighbours" << std::endl;

               for(int i = 0; i < num_materials; i++){
                  exchange_matrix_array.resize(num_materials);
                  for(int j = 0; j < num_materials; j++){
                     exchange_matrix_array[i].resize(num_materials);
                     for(int k = 0; k < num_neighbours; k++) exchange_matrix_array[i][j].resize(num_neighbours);
                  }
               }
               // save maximum array sizes
               max_materials  = num_materials;
               max_neighbours = num_neighbours;
               return;
            }

         public:

            // constructor
            exchange_matrix_4D_t(){
               // set initial sizes of zero
               max_materials = 0;
               max_neighbours = 0;
               // resize storage
               resize(max_materials, max_neighbours);
            };

            // function to set a particular set of exchange values
            void set_exchange_values(int material_i, int material_j, int neighbour, std::vector<double> exchange_values){

               // check for expansion of storage size
               int new_max_materials = max_materials;
               if( material_i > new_max_materials-1 ) new_max_materials = material_i+1;
               if( material_j > new_max_materials-1 ) new_max_materials = material_j+1;

               int new_max_neighbours = max_neighbours;
               if( neighbour > new_max_neighbours-1 ) new_max_neighbours = neighbour+1;

               if( new_max_materials > max_materials || new_max_neighbours > max_neighbours) resize(new_max_materials, new_max_neighbours);

               // now add exchange values in-place
               exchange_matrix_array[material_i][material_j][neighbour] = exchange_values;

               return;

            }

            // function to get a list of exchange values for a certain material pair and neighbour number
            std::vector<double> get_exchange_values(int material_i, int material_j, int neighbour){

               // initialise a vector to return exchange values
               std::vector<double> exchange_values(1,0.0);

               // check for valid array access
               bool ok = true;
               if(material_i > max_materials-1) ok = false;
               if(material_j > max_materials-1) ok = false;
               if(neighbour > max_neighbours-1) ok = false;

               // if out of bounds access, this means constant was never set so assume zero
               if(!ok) return exchange_values;

               // array index value is valid, so return values (single value, Jij=0.0 array if empty)
               if(exchange_matrix_array[material_i][material_j][neighbour].size() > 0){
                  exchange_values = exchange_matrix_array[material_i][material_j][neighbour];
               }

               // return exchange values
               return exchange_values;

            }

            int get_max_shells(){
               return max_neighbours;
            }

            int get_max_materials(){
               return max_materials;
            }

      };

      //-----------------------------------------------------------------------------
      // materials class for storing exchange material parameters
      //-----------------------------------------------------------------------------
      class mp_t{

         private:

         public:

            // variables
            std::vector<double> dmi; // Dzyaloshinskii-Moriya interaction constant
            std::vector<double> kitaev; // Dzyaloshinskii-Moriya interaction constant

            // constructor
            mp_t (const unsigned int max_materials = 100)
            {
               // resize arrays to correct size
               dmi.resize(max_materials, 0.0); // initialise pair anisotropy constants to zero
               kitaev.resize(max_materials, 0.0); // initialise pair anisotropy constants to zero

            }; // end of constructor

      }; // end of exchange::internal::mp class

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------
      extern std::vector<internal::mp_t> mp; // array of material properties

      extern exchange_matrix_4D_t bilinear_exchange_constants; // array of exchange constants
      extern exchange_matrix_4D_t biquadratic_exchange_constants; // array of biquadratic exchange constants

      extern bool enable_dmi; // flag to enable dmi calculation
      extern bool enable_kitaev; // flag to enable kitaev calculation

      extern double dmi_cutoff_range;    // cutoff range for DMI calculation (Ångstroms)
      extern double kitaev_cutoff_range; // cutoff range for Kitaev calculation (Ångstroms)
      extern double exchange_factor;     // scaling factor for exchange constants (usually to correct for ab-initio)

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
      void calculate_kitaev(std::vector<std::vector <neighbours::neighbour_t> >& cneighbourlist);
      void unroll_exchange_interactions(std::vector<std::vector <neighbours::neighbour_t> >& bilinear);
      void unroll_normalised_exchange_interactions(std::vector<std::vector <neighbours::neighbour_t> >& bilinear);
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
