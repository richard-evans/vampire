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
#include "cusp/print.h"


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

            spin3N.assign( 3*::atoms::num_atoms, 0);
            field3N.assign( 3*::atoms::num_atoms, 0);

            //Local storage for nbr list
            std::vector<int> row_inds;
            std::vector<int> col_inds;
            std::vector<double> vals;

            // tolerance to ignore exchange components
            const double tol = 1e-10;
            // loop over all atoms
            for(int atom = 0; atom < ::atoms::num_atoms; ++atom){

                // temporray constants for loop start and end indices
                const int start = ::atoms::neighbour_list_start_index[atom];
                const int end   = ::atoms::neighbour_list_end_index[atom]+1;

                // loop over all neighbours
                for(int nn = start; nn < end; ++nn){

                    const int natom = ::atoms::neighbour_list_array[nn]; // get neighbouring atom number
                    const int iid = ::atoms::neighbour_interaction_type_array[nn]; // interaction id

                    // Store nbr list in local array
                    // Nbr list is expanded to 3N by 3N size
                    switch( ::exchange::get_exchange_type())
                    {
                        case 0: // Isotropic
                            for ( int a = 0; a < 3; a++){
                                double Jab = ::atoms::i_exchange_list[iid].Jij;
                                if( fabs(Jab) > tol) {
                                    row_inds.push_back( atom + a*Natoms);
                                    col_inds.push_back( natom + a*Natoms);
                                    vals.push_back( Jab);
                                }
                            }
                            break;

                        case 1: // Vectorial
                            for ( int a = 0; a < 3; a++){
                                double Jab = ::atoms::v_exchange_list[iid].Jij[a];
                                if( fabs(Jab) > tol) {
                                    row_inds.push_back( atom + a*Natoms);
                                    col_inds.push_back( natom + a*Natoms);
                                    vals.push_back( Jab);
                                }
                            }
                            break;

                        case 2: // Tensor
                            for ( int a = 0; a < 3; a++){
                                for ( int b = 0; b < 3; b++){
                                    double Jab = ::atoms::t_exchange_list[iid].Jij[a][b];
                                    if( fabs(Jab) > tol) {
                                        row_inds.push_back( atom + a*Natoms);
                                        col_inds.push_back( natom + b*Natoms);
                                        vals.push_back( Jab);
                                    }
                                }
                            }
                            break;
                    }

                }

            }


            // Declare a local matrix on the host using coordinate format to be filled
            cusp::coo_matrix< int, cu::cu_real_t, cusp::host_memory> J_matrix_h;


            zlog << zTs() << "Expanded CPU nbr list to 3N by 3N format, no. of non-zeros is :" << vals.size() << " with a tol = " << tol << std::endl;

            // Set COO matrix to size 3Natoms by 3Natoms with number of non-zeros found
            J_matrix_h.resize(
                    3*::atoms::num_atoms,
                    3*::atoms::num_atoms,
                    vals.size()
                    );

            // copy in to CUSP COO matrix for easy conversion
            for( int i = 0; i < vals.size(); i++) {
                J_matrix_h.row_indices[i] = row_inds[i];
                J_matrix_h.column_indices[i] = col_inds[i];
                J_matrix_h.values[i] = vals[i];
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
            switch( ::exchange::get_exchange_type())
            {
                case 0: // Isotropic
                    J_isot_initialised = true;
                    break;
                case 1: // Vectorial
                    J_vect_initialised = true;
                    break;
                case 2: // Tensor
                    J_tens_initialised = true;
                    break;
            }

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

            check_cuda_errors(__FILE__,__LINE__);
            cusp::multiply(
                  J_matrix_d,
                  spin3N,
                  field3N);

            check_cuda_errors(__FILE__,__LINE__);
            thrust::copy( field3N.begin(), field3N.begin() + ::atoms::num_atoms, cu::x_total_spin_field_array.begin() );
            thrust::copy( field3N.begin() + ::atoms::num_atoms, field3N.begin() + 2*::atoms::num_atoms, cu::y_total_spin_field_array.begin() );
            thrust::copy( field3N.begin() + 2*::atoms::num_atoms, field3N.end(), cu::z_total_spin_field_array.begin() );


            check_cuda_errors(__FILE__,__LINE__);
            return EXIT_SUCCESS;
         }
      } // end namespace exchange

   } // end namespace internal

} // end namespace vcuda
