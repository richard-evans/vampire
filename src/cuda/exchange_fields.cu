#include "exchange_fields.hpp"

#include "atoms.hpp"

#include "cuda_utils.hpp"
#include "internal.hpp"
#include "data.hpp"

#include <vector>


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


         int initialise_exchange()
         {

            check_device_memory(__FILE__,__LINE__);

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
                  check_cuda_errors(__FILE__,__LINE__);

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
                  check_cuda_errors(__FILE__,__LINE__);

                  break;

               case 2: // Tensor
                  break;
            }

            exchange_initialised = true;

            std::cout << "Made matrix" << std::endl;
            check_device_memory(__FILE__,__LINE__);
            check_cuda_errors(__FILE__,__LINE__);
            return EXIT_SUCCESS;
         }


         int finalise_exchange()
         {
            check_cuda_errors(__FILE__,__LINE__);
            return EXIT_SUCCESS;
         }

         int calculate_exchange_fields()
         {

            check_cuda_errors(__FILE__,__LINE__);

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

