#include "exchange_fields.hpp"

#include "atoms.hpp"
#include "exchange.hpp"
#include "vio.hpp"
#include "cuda_utils.hpp"
#include "internal.hpp"
#include "data.hpp"

#include <vector>


#include "cusp/array2d.h"
#include "cusp/coo_matrix.h"


int calculate_exchange_fields(int, int);

namespace cu = vcuda::internal;

namespace vcuda
{

   namespace internal
   {

      namespace exchange
      {


         bool exchange_initialised = false;

         bool J_isot_initialised = false;
         bool J_vect_initialised = false;
         bool J_tens_initialised = false;

         cu_real_array_t   spin3N;
         cu_real_array_t   field3N;

         cu_exch_mat_t  J_matrix_d;

         int initialise_exchange()
         {

            // print out informative message regarding compile time option for matrix format
            #if CUDA_MATRIX == CSR
               zlog << zTs() << "Configured exchange calculation using CSR matrix format" << std::endl;
            #elif CUDA_MATRIX == DIA
               zlog << zTs() << "Configured exchange calculation using DIA matrix format" << std::endl;
            #elif CUDA_MATRIX == ELL
               zlog << zTs() << "Configured exchange calculation using ELL matrix format" << std::endl;
            #else
               zlog << zTs() << "Configured exchange calculation using default CSR matrix format" << std::endl;
            #endif

            check_device_memory(__FILE__,__LINE__);

            const int Natoms = ::atoms::num_atoms;
            const int Nnbrs = ::atoms::neighbour_list_array.size();

            spin3N.assign( 3*::atoms::num_atoms, 0);
            field3N.assign( 3*::atoms::num_atoms, 0);

            cusp::array1d<int, cusp::host_memory> row_indices( Nnbrs );
            cusp::array1d<int, cusp::host_memory> row_offsets( Natoms + 1);
            row_offsets[0] = 0;
            for( int atom = 0; atom < Natoms; atom++) row_offsets[atom+1] = ::atoms::neighbour_list_end_index[atom]+1;
            cusp::offsets_to_indices( row_offsets, row_indices);

            // Declare a local matrix on the host using coordinate format to be filled
            cusp::coo_matrix< int, cu::cu_real_t, cusp::host_memory> J_matrix_h;

            switch( ::atoms::exchange_type)
            {
                case 0: // Isotropic
                    J_matrix_h.resize(
                            3*::atoms::num_atoms,
                            3*::atoms::num_atoms,
                            3*::atoms::neighbour_list_array.size()
                            );
                    break;

                case 1: // Vectorial
                    J_matrix_h.resize(
                            3*::atoms::num_atoms,
                            3*::atoms::num_atoms,
                            3*::atoms::neighbour_list_array.size()
                            );
                    break;

                case 2: // Tensor
                    J_matrix_h.resize(
                            9*::atoms::num_atoms,
                            9*::atoms::num_atoms,
                            9*::atoms::neighbour_list_array.size()
                            );
                    break;
            }



            for( int i = 0; i < Nnbrs; i++)
            {


               int iid = ::atoms::neighbour_interaction_type_array[i];

               switch(exchange_type)
               {
                   case 0: // Isotropic
                       J_matrix_h.row_indices[i]            = row_indices[i];
                       J_matrix_h.column_indices[i]         = ::atoms::neighbour_list_array[i];
                       J_matrix_h.values[i]                 = - ::atoms::i_exchange_list[iid].Jij;

                       J_matrix_h.row_indices[i+Nnbrs]      = row_indices[i]+Natoms;
                       J_matrix_h.column_indices[i+Nnbrs]   = ::atoms::neighbour_list_array[i] + Natoms;
                       J_matrix_h.values[i+Nnbrs]           = - ::atoms::i_exchange_list[iid].Jij;

                       J_matrix_h.row_indices[i+2*Nnbrs]    = row_indices[i]+2*Natoms;
                       J_matrix_h.column_indices[i+2*Nnbrs] = ::atoms::neighbour_list_array[i] + 2*Natoms;
                       J_matrix_h.values[i+2*Nnbrs]         = - ::atoms::i_exchange_list[iid].Jij;
                       break;

                   case 1: // Vectorial
                       J_matrix_h.row_indices[i]            = row_indices[i];
                       J_matrix_h.column_indices[i]         = ::atoms::neighbour_list_array[i];
                       J_matrix_h.values[i]                 = - ::atoms::v_exchange_list[iid].Jij[0];

                       J_matrix_h.row_indices[i+Nnbrs]      = row_indices[i]+Natoms;
                       J_matrix_h.column_indices[i+Nnbrs]   = ::atoms::neighbour_list_array[i] + Natoms;
                       J_matrix_h.values[i+Nnbrs]           = - ::atoms::v_exchange_list[iid].Jij[1];

                       J_matrix_h.row_indices[i+2*Nnbrs]    = row_indices[i]+2*Natoms;
                       J_matrix_h.column_indices[i+2*Nnbrs] = ::atoms::neighbour_list_array[i] + 2*Natoms;
                       J_matrix_h.values[i+2*Nnbrs]         = - ::atoms::v_exchange_list[iid].Jij[2];
                       break;

                   case 2: // Tensor

                       J_matrix_h.row_indices[i]            = row_indices[i];
                       J_matrix_h.column_indices[i]         = ::atoms::neighbour_list_array[i];
                       J_matrix_h.values[i]                 = - ::atoms::t_exchange_list[iid].Jij[0][0];

                       J_matrix_h.row_indices[i+1]            = row_indices[i];
                       J_matrix_h.column_indices[i+1]         = ::atoms::neighbour_list_array[i]+Natoms;
                       J_matrix_h.values[i+1]                 = - ::atoms::t_exchange_list[iid].Jij[0][1];

                       J_matrix_h.row_indices[i+2]            = row_indices[i];
                       J_matrix_h.column_indices[i+2]         = ::atoms::neighbour_list_array[i]+2*Natoms;
                       J_matrix_h.values[i+2]                 = - ::atoms::t_exchange_list[iid].Jij[0][2];

                       J_matrix_h.row_indices[i+Nnbrs]      = row_indices[i]+Natoms;
                       J_matrix_h.column_indices[i+Nnbrs]   = ::atoms::neighbour_list_array[i];
                       J_matrix_h.values[i+Nnbrs]           = - ::atoms::t_exchange_list[iid].Jij[1][0];

                       J_matrix_h.row_indices[i+Nnbrs+1]      = row_indices[i]+Natoms;
                       J_matrix_h.column_indices[i+Nnbrs+1]   = ::atoms::neighbour_list_array[i] + Natoms;
                       J_matrix_h.values[i+Nnbrs+1]           = - ::atoms::t_exchange_list[iid].Jij[1][1];

                       J_matrix_h.row_indices[i+Nnbrs+2]      = row_indices[i]+Natoms;
                       J_matrix_h.column_indices[i+Nnbrs+2]   = ::atoms::neighbour_list_array[i] + 2*Natoms;
                       J_matrix_h.values[i+Nnbrs+2]           = - ::atoms::t_exchange_list[iid].Jij[1][2];

                       J_matrix_h.row_indices[i+2*Nnbrs]    = row_indices[i]+2*Natoms;
                       J_matrix_h.column_indices[i+2*Nnbrs] = ::atoms::neighbour_list_array[i];
                       J_matrix_h.values[i+2*Nnbrs]         = - ::atoms::t_exchange_list[iid].Jij[2][0];

                       J_matrix_h.row_indices[i+2*Nnbrs+1]    = row_indices[i]+2*Natoms;
                       J_matrix_h.column_indices[i+2*Nnbrs+1] = ::atoms::neighbour_list_array[i] + Natoms;
                       J_matrix_h.values[i+2*Nnbrs+1]         = - ::atoms::t_exchange_list[iid].Jij[2][1];

                       J_matrix_h.row_indices[i+2*Nnbrs+2]    = row_indices[i]+2*Natoms;
                       J_matrix_h.column_indices[i+2*Nnbrs+2] = ::atoms::neighbour_list_array[i] + 2*Natoms;
                       J_matrix_h.values[i+2*Nnbrs+2]         = - ::atoms::t_exchange_list[iid].Jij[2][2];

                       break;
               }
            }

            // Use the sort member function to double check ordering before convert
            J_matrix_h.sort_by_row_and_column();

            // Black magic to turn CUDA_MATRIX into a string
            #define STRING(s) #s
            #define VALSTRING(s) STRING(s)

            // Print out informative message
            zlog << zTs() << "Attempting matrix conversion from COO to " << VALSTRING(CUDA_MATRIX) << " and transferring to device..." << std::endl;

            cusp::convert( J_matrix_h, J_matrix_d);

            zlog << zTs() << "Matrix conversion and transfer complete." << std::endl;

            exchange_initialised = true;

            check_device_memory(__FILE__,__LINE__);
            check_cuda_errors(__FILE__,__LINE__);
            return EXIT_SUCCESS;
         }


         int finalise_exchange()
         {
            spin3N.cu_real_array_t::~cu_real_array_t();
            field3N.cu_real_array_t::~cu_real_array_t();
            J_matrix_d.cu_exch_mat_t::~cu_exch_mat_t ();

            check_cuda_errors(__FILE__,__LINE__);
            return EXIT_SUCCESS;
         }

         int calculate_exchange_fields()
         {

            check_cuda_errors(__FILE__,__LINE__);

            if( !exchange_initialised) initialise_exchange();

            thrust::copy( cu::atoms::x_spin_array.begin(), cu::atoms::x_spin_array.end(), spin3N.begin());
            thrust::copy( cu::atoms::y_spin_array.begin(), cu::atoms::y_spin_array.end(), spin3N.begin() + ::atoms::num_atoms);
            thrust::copy( cu::atoms::z_spin_array.begin(), cu::atoms::z_spin_array.end(), spin3N.begin() + 2*::atoms::num_atoms);

            cusp::multiply(
                  J_matrix_d,
                  spin3N,
                  field3N);

            thrust::copy( field3N.begin(), field3N.begin() + ::atoms::num_atoms, cu::x_total_spin_field_array.begin() );
            thrust::copy( field3N.begin() + ::atoms::num_atoms, field3N.begin() + 2*::atoms::num_atoms, cu::y_total_spin_field_array.begin() );
            thrust::copy( field3N.begin() + 2*::atoms::num_atoms, field3N.end(), cu::z_total_spin_field_array.begin() );


            check_cuda_errors(__FILE__,__LINE__);
            return EXIT_SUCCESS;
         }
      } // end namespace exchange

   } // end namespace internal

} // end namespace vcuda
