#include "exchange_fields.hpp"

#include "atoms.hpp"
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

            spin3N.assign( 3*::atoms::num_atoms, 0);
            field3N.assign( 3*::atoms::num_atoms, 0);

            cusp::csr_matrix < int, cu::cu_real_t, cusp::host_memory > J_xx_matrix_h (
                  ::atoms::num_atoms,
                  ::atoms::num_atoms,
                  ::atoms::neighbour_list_array.size()
                  );

            cusp::csr_matrix < int, cu::cu_real_t, cusp::host_memory > J_yy_matrix_h (
                  ::atoms::num_atoms,
                  ::atoms::num_atoms,
                  ::atoms::neighbour_list_array.size()
                  );

            cusp::csr_matrix < int, cu::cu_real_t, cusp::host_memory > J_zz_matrix_h (
                  ::atoms::num_atoms,
                  ::atoms::num_atoms,
                  ::atoms::neighbour_list_array.size()
                  );


            J_xx_matrix_h.row_offsets[0] = 0.0;
            J_yy_matrix_h.row_offsets[0] = 0.0;
            J_zz_matrix_h.row_offsets[0] = 0.0;

            for (int atom = 0; atom < ::atoms::num_atoms; atom++) {
               J_xx_matrix_h.row_offsets[atom+1] = ::atoms::neighbour_list_end_index[atom]+1;
               J_yy_matrix_h.row_offsets[atom+1] = ::atoms::neighbour_list_end_index[atom]+1;
               J_zz_matrix_h.row_offsets[atom+1] = ::atoms::neighbour_list_end_index[atom]+1;
            }

            for (int i = 0; i < ::atoms::neighbour_list_array.size(); i++) {
               J_xx_matrix_h.column_indices[i] = ::atoms::neighbour_list_array[i];
               J_yy_matrix_h.column_indices[i] = ::atoms::neighbour_list_array[i];
               J_zz_matrix_h.column_indices[i] = ::atoms::neighbour_list_array[i];
            }

            cusp::array1d<int, cusp::host_memory> row_indices( J_xx_matrix_h.column_indices);
            cusp::array1d<int, cusp::host_memory> column_indices( J_xx_matrix_h.column_indices);
            cusp::offsets_to_indices( J_xx_matrix_h.row_offsets, row_indices);

            cusp::coo_matrix< int, cu::cu_real_t, cusp::host_memory> J_matrix_h(
                  3*::atoms::num_atoms,
                  3*::atoms::num_atoms,
                  3*::atoms::neighbour_list_array.size()
                  );

            //std::cerr << "Attempting to fill matrix." << std::endl;

            const int Natoms = ::atoms::num_atoms;
            const int Nnbrs = row_indices.size();
            for( int i = 0; i < Nnbrs; i++)
            {
               J_matrix_h.row_indices[i] = row_indices[i];
               J_matrix_h.row_indices[i+Nnbrs] = row_indices[i]+Natoms;
               J_matrix_h.row_indices[i+2*Nnbrs] = row_indices[i]+2*Natoms;

               J_matrix_h.column_indices[i] = column_indices[i];
               J_matrix_h.column_indices[i+Nnbrs] = column_indices[i]+Natoms;
               J_matrix_h.column_indices[i+2*Nnbrs] = column_indices[i]+2*Natoms;

               int iid = ::atoms::neighbour_interaction_type_array[i];
               switch( ::atoms::exchange_type)
               {
                  case 0: // Isotropic
                     J_matrix_h.values[i]         = - ::atoms::i_exchange_list[iid].Jij;
                     J_matrix_h.values[i+Nnbrs]   = - ::atoms::i_exchange_list[iid].Jij;
                     J_matrix_h.values[i+2*Nnbrs] = - ::atoms::i_exchange_list[iid].Jij;
                     break;

                  case 1: // Vectorial
                     J_matrix_h.values[i]         = - ::atoms::v_exchange_list[iid].Jij[0];
                     J_matrix_h.values[i+Nnbrs]   = - ::atoms::v_exchange_list[iid].Jij[1];
                     J_matrix_h.values[i+2*Nnbrs] = - ::atoms::v_exchange_list[iid].Jij[2];
                     break;

                  case 2: // Tensor
                     std::cerr << "Error! Tensorial form of exchange not yet implemented in cuda version!" << std::endl;
                     zlog << zTs() << "Error! Tensorial form of exchange not yet implemented in cuda version!" << std::endl;
                     break;
               }
            }

            // Black magic to turn CUDA_MATRIX into a string
            #define STRING(s) #s
            #define VALSTRING(s) STRING(s)

            // Print out informative message
            zlog << zTs() << "Attempting matrix conversion from CSR to " << VALSTRING(CUDA_MATRIX) << " and transferring to device..." << std::endl;

            cusp::convert( J_matrix_h, J_matrix_d);

            zlog << zTs() << "Matrix conversion and transfer complete." << std::endl;

            const size_t occupied_diagonals = count_diagonals(J_xx_matrix_h.num_rows, J_xx_matrix_h.num_rows, row_indices, J_xx_matrix_h.column_indices);
            const float size       = float(occupied_diagonals) * float(J_xx_matrix_h.num_rows);
            const float fill_ratio = size / std::max(1.0f, float(J_xx_matrix_h.num_entries));

            zlog << zTs() << "Cuda Matrix: Diagonals = " << occupied_diagonals << ", size = " << size << ", fill ratio = "<< fill_ratio << std::endl;

            switch( ::atoms::exchange_type)
            {
               case 0: // Isotropic

                  //--------------------------------------------------------------
                  // Exchange is isotropic so Jxx = Jyy = Jzz
                  // and Jxy = Jxz = Jyx = 0
                  //--------------------------------------------------------------

                  for (int i = 0; i < ::atoms::neighbour_list_array.size(); i++) {
                     int iid = ::atoms::neighbour_interaction_type_array[i];
                     J_xx_matrix_h.values[i] = - ::atoms::i_exchange_list[iid].Jij;
                  }



                  cusp::convert(J_xx_matrix_h, J_xx_mat_d);

                  J_isot_initialised = true;

                  break;

               case 1: // Vector

                  //--------------------------------------------------------------
                  // Exchange is diagonal so Jxx != Jyy != Jzz
                  // and Jxy = Jxz = Jyx = 0
                  //--------------------------------------------------------------

                  // Copy J values from vampire exchange list to values list
                  for( int i = 0; i < ::atoms::neighbour_list_array.size(); i++)
                  {
                     int iid = ::atoms::neighbour_interaction_type_array[i];
                     J_xx_matrix_h.values[i] = - ::atoms::v_exchange_list[iid].Jij[0];
                     J_yy_matrix_h.values[i] = - ::atoms::v_exchange_list[iid].Jij[1];
                     J_zz_matrix_h.values[i] = - ::atoms::v_exchange_list[iid].Jij[2];
                  }

                  cusp::convert(J_xx_matrix_h, J_xx_mat_d);
                  cusp::convert(J_yy_matrix_h, J_yy_mat_d);
                  cusp::convert(J_zz_matrix_h, J_zz_mat_d);

                  J_vect_initialised = true;

                  break;

               case 2: // Tensor
                  std::cerr << "Error! Tensorial form of exchange not yet implemented in cuda version!" << std::endl;
                  zlog << zTs() << "Error! Tensorial form of exchange not yet implemented in cuda version!" << std::endl;
                  break;
            }

            exchange_initialised = true;

            //std::cout << "Made matrix" << std::endl;
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
              /*
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

              */
            switch( ::atoms::exchange_type)
            {
               case 0: // Isotropic

                  //--------------------------------------------------------------
                  // Exchange is isotropic so Jxx = Jyy = Jzz
                  // and Jxy = Jxz = Jyx = 0
                  //--------------------------------------------------------------

                  if( !exchange_initialised) initialise_exchange();

                  // FIXME This maybe boosted
                  // It should keep the old values stored in the spin field
                  // Since Jxx = Jyy = Jzz only the Jxx array is used

                  cusp::multiply(
                        J_xx_mat_d,
                        cu::atoms::x_spin_array,
                        cu::x_total_spin_field_array
                        );
                  cusp::multiply(
                        J_xx_mat_d,
                        cu::atoms::y_spin_array,
                        cu::y_total_spin_field_array
                        );
                  cusp::multiply(
                        J_xx_mat_d,
                        cu::atoms::z_spin_array,
                        cu::z_total_spin_field_array
                        );

                  break;

               case 1: // Vector exchange

                  //--------------------------------------------------------------
                  // Exchange is diagonal so Jxx != Jyy != Jzz
                  // and Jxy = Jxz = Jyx = 0
                  //--------------------------------------------------------------

                  if( !exchange_initialised) initialise_exchange();

                  cusp::multiply(
                        J_xx_mat_d,
                        cu::atoms::x_spin_array,
                        cu::x_total_spin_field_array
                        );
                  cusp::multiply(
                        J_yy_mat_d,
                        cu::atoms::y_spin_array,
                        cu::y_total_spin_field_array
                        );
                  cusp::multiply(
                        J_zz_mat_d,
                        cu::atoms::z_spin_array,
                        cu::z_total_spin_field_array
                        );

                  break;
            }

            check_cuda_errors(__FILE__,__LINE__);
            return EXIT_SUCCESS;
         }
      } // end namespace exchange

   } // end namespace internal

} // end namespace vcuda
