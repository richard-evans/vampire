#include "exchange_fields.hpp"

#include "atoms.hpp"
#include "exchange.hpp"
#include "vio.hpp"
#include "cuda_utils.hpp"
#include "internal.hpp"
#include "data.hpp"

#include <vector>

#include "cusparse.h"

int calculate_exchange_fields(int, int);

namespace cu = vcuda::internal;

namespace vcuda
{

   namespace internal
   {

      namespace exchange
      {


         bool exchange_initialised = false;
         bool empty_exchange = true;

         bool J_isot_initialised = false;
         bool J_vect_initialised = false;
         bool J_tens_initialised = false;

        cu_real_t *d_spin3n;
        cu_real_t *d_field3n;

        //cu_real_array_t   spin3N;
        //cu_real_array_t   field3N;

         cusparseDnVecDescr_t vecX, vecY;
         cu_exch_mat_t  J_matrix_d;

         // Arrays to store the coo and csr matrix
        int *d_csr_rows;
        int *d_coo_rows;
        int *d_coo_cols;

        cu_real_t *d_coo_vals;

        //cu_index_array_t csr_rows_d;
        //cu_index_array_t coo_rows_d;
        //cu_index_array_t coo_cols_d;
        //cu_real_array_t  coo_vals_d;

         size_t Nrows = 0;
         size_t Ncols = 0;
         size_t Nnz = 0;

         // CUSPARSE contex handle and error
         cusparseHandle_t handle=0;
         cusparseStatus_t status;

         // Constants for the matrix-vector product y = alpha*A*x + beta*y
         cu_real_t alpha = 1.0;
         cu_real_t beta = 0.0;

         // Buffer workspace for the spmv product
         void * spmv_buffer_d = NULL;
         size_t buffer_size = 0;

         // Routine to sort the coo sparse matrix by the row major order index
         void sort_coo_list(std::vector<int> &rows, std::vector<int> &cols, std::vector<cu_real_t> &vals, const size_t Nrows, const size_t Ncols )
         {
             // Create a list of id-value pairs
             std::vector< std::pair<size_t, cu_real_t> > list;

             // Calculate the row major order index and create pair list
             for ( size_t i = 0; i < vals.size(); i++ ) {
                 size_t id = cols[i] + rows[i]*Ncols;
                 list.push_back( std::make_pair(id, vals[i]));
             }

             std::sort( list.begin(), list.end() );

             // Overwrite the input vectors with sorted list
             for( size_t i = 0; i < list.size(); i++) {
                 size_t id = list[i].first;
                 size_t row = size_t( id / Ncols);
                 size_t col = id % Ncols;
                 // Is that ok to do, or are we going to truncate?
                 rows[i] = row;
                 cols[i] = col;
                 vals[i] = list[i].second;
             }


         }

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

            //cudaMalloc((void**)&d_spin3n, 3 * ::atoms::num_atoms * sizeof(cu_real_t));
            //cudaMalloc((void**)&d_field3n, 3 * ::atoms::num_atoms * sizeof(cu_real_t));

            // NOTE: This will NOT work for values other than 0
            //cudaMemset(d_spin3n, 0, 3 * ::atoms::num_atoms * sizeof(cu_real_t));
            //cudaMemset(d_field3n, 0, 3 * ::atoms::num_atoms * sizeof(cu_real_t));

            //spin3N.assign( 3*::atoms::num_atoms, 0);
            //field3N.assign( 3*::atoms::num_atoms, 0);

            //Local storage for nbr list
            std::vector<int> row_inds;
            std::vector<int> col_inds;
            std::vector<cu_real_t> vals;

            // NOTE: Below calculations have been left in double precision
            // for now !!
            // tolerance to ignore exchange components
            const cu_real_t tol = 1e-10;
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
                                cu_real_t Jab = ::atoms::i_exchange_list[iid].Jij;
                                if( fabs(Jab) > tol) {
                                    row_inds.push_back( atom + a*Natoms);
                                    col_inds.push_back( natom + a*Natoms);
                                    vals.push_back( Jab);
                                }
                            }
                            break;

                        case 1: // Vectorial
                            for ( int a = 0; a < 3; a++){
                                cu_real_t Jab = ::atoms::v_exchange_list[iid].Jij[a];
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
                                    cu_real_t Jab = ::atoms::t_exchange_list[iid].Jij[a][b];
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

            zlog << zTs() << "Expanded CPU nbr list to 3N by 3N format, no. of non-zeros is :" << vals.size() << " with a tol = " << tol << std::endl;

            Ncols = 3 * Natoms;
            Nrows = 3 * Natoms;
            Nnz = vals.size();

            if ( Nnz > 0) {


                // sort coo list to make sure the data is sequential
                sort_coo_list(row_inds, col_inds, vals, Nrows, Ncols);

                // allocate space for the device data
                /*
                coo_rows_d.resize(Nnz);
                coo_cols_d.resize(Nnz);
                coo_vals_d.resize(Nnz);
                csr_rows_d.resize(Nrows+1);
                */
                cudaMalloc((void**)&d_coo_rows, Nnz * sizeof(int));
                cudaMalloc((void**)&d_coo_cols, Nnz * sizeof(int));
                cudaMalloc((void**)&d_csr_rows, (Nrows + 1) * sizeof(int));
                cudaMalloc((void**)&d_coo_vals, Nnz * sizeof(cu_real_t));

                cudaMemcpy(d_coo_rows, row_inds.data(), Nnz * sizeof(int), cudaMemcpyHostToDevice);
                cudaMemcpy(d_coo_cols, col_inds.data(), Nnz * sizeof(int), cudaMemcpyHostToDevice);
                cudaMemcpy(d_coo_vals, vals.data(), Nnz * sizeof(cu_real_t), cudaMemcpyHostToDevice);

                // initialise cusparse handle
                status = cusparseCreate(&handle);
                if (status != CUSPARSE_STATUS_SUCCESS) {
                    std::cerr << "Failed to initialise CUSPARSE library!" << std::endl;
                    return 1;
                }

                // cusparse routine to convert coo row data into csr row offsets
                status = cusparseXcoo2csr(  handle,
                                            d_coo_rows,
                                            Nnz,
                                            Ncols,
                                            d_csr_rows,
                                            CUSPARSE_INDEX_BASE_ZERO);

                // create the CSR descriptor for CUSPARSE
                status = cusparseCreateCsr(&J_matrix_d, Nrows, Ncols, Nnz,
                                            d_csr_rows,
                                            d_coo_cols,
                                            d_coo_vals,
                                          CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                          CUSPARSE_INDEX_BASE_ZERO, CUSPARSE_REAL);

                 if (status != CUSPARSE_STATUS_SUCCESS) {
                     std::cerr << "Failed to initialise sparse matrix descriptor!" << std::endl;
                     return 1;
                 }



                 std::vector<int> csr_row_inds;
                 csr_row_inds.resize(Nrows+1);
                 cudaMemcpy(csr_row_inds.data(), d_csr_rows, (Nrows+1) * sizeof(int), cudaMemcpyDeviceToHost);

                 //for (int i = 0; i < Nrows; i++){
                 //    for (int j = csr_row_inds[i]; j < csr_row_inds[i+1]; j++){
                 //        std::cerr << i << "  " << col_inds[j] << "  " << vals[j] << std::endl;
                 //    }
                 //}

                 // Create the dense vector descriptors for the input and output (Y = A*X)
                 //cusparseCreateDnVec( &vecX, Ncols, d_spin3n, CUSPARSE_REAL );
                 //cusparseCreateDnVec( &vecY, Nrows, d_field3n, CUSPARSE_REAL );
                 cusparseCreateDnVec( &vecX, Ncols, cu::atoms::d_spin, CUSPARSE_REAL );
                 cusparseCreateDnVec( &vecY, Nrows, cu::d_spin_field, CUSPARSE_REAL );

                // allocate an external buffer if needed
                status = cusparseSpMV_bufferSize(handle,
                                    CUSPARSE_OPERATION_NON_TRANSPOSE,
                                    &alpha, J_matrix_d, vecX, &beta, vecY,
                                    CUSPARSE_REAL,
                                    CUSPARSE_CSRMV_ALG1,
                                    &buffer_size);

                cudaMalloc(&spmv_buffer_d, buffer_size);

                // Declare a local matrix on the host using coordinate format to be filled
                //cusp::coo_matrix< int, cu::cu_real_t, cusp::host_memory> J_matrix_h;

                // Set COO matrix to size 3Natoms by 3Natoms with number of non-zeros found
                //J_matrix_h.resize(
                //        3*::atoms::num_atoms,
                //        3*::atoms::num_atoms,
                //        vals.size()
                //        );

                // copy in to CUSP COO matrix for easy conversion
                //for( int i = 0; i < vals.size(); i++) {
                //    J_matrix_h.row_indices[i] = row_inds[i];
                //    J_matrix_h.column_indices[i] = col_inds[i];
                //    J_matrix_h.values[i] = vals[i];
                //}

                // Use the sort member function to double check ordering before convert
                //J_matrix_h.sort_by_row_and_column();

                // Black magic to turn CUDA_MATRIX into a string
                #define STRING(s) #s
                #define VALSTRING(s) STRING(s)

                // Print out informative message
                //zlog << zTs() << "Attempting matrix conversion from COO to " << VALSTRING(CUDA_MATRIX) << " and transferring to device..." << std::endl;

                //cusp::convert( J_matrix_h, J_matrix_d);

                zlog << zTs() << "Matrix conversion and transfer complete." << std::endl;

                empty_exchange = false;
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
            } // End if number of non-zeros greater than zero

            exchange_initialised = true;

            check_device_memory(__FILE__,__LINE__);
            check_cuda_errors(__FILE__,__LINE__);
            return EXIT_SUCCESS;
         }


         int finalise_exchange()
         {
            //cudaFree(d_spin3n);
            //cudaFree(d_field3n);
            //spin3N.cu_real_array_t::~cu_real_array_t();
            //field3N.cu_real_array_t::~cu_real_array_t();
            //J_matrix_d.cu_exch_mat_t::~cu_exch_mat_t ();

            cudaFree(d_csr_rows);
            cudaFree(d_coo_rows);
            cudaFree(d_coo_cols);
            cudaFree(d_coo_vals);

            //csr_rows_d.cu_index_array_t::~cu_index_array_t();
            //coo_rows_d.cu_index_array_t::~cu_index_array_t();
            //coo_cols_d.cu_index_array_t::~cu_index_array_t();
            //coo_vals_d.cu_real_array_t::~cu_real_array_t();

            if( !empty_exchange ) {
                // destroy vector and matrix descriptors
                cusparseDestroyDnVec(vecX);
                cusparseDestroyDnVec(vecY);
                status = cusparseDestroySpMat(J_matrix_d);

                if (status != CUSPARSE_STATUS_SUCCESS) {
                    std::cerr << "Matrix descriptor destruction failed" << std::endl;
                    return 1;
                }

                // destroy cusparse handle
                status = cusparseDestroy(handle);
                handle = 0;
                if (status != CUSPARSE_STATUS_SUCCESS) {
                    std::cerr << "CUSPARSE Library release of resources failed" << std::endl;
                    return 1;
                }
            }

            check_cuda_errors(__FILE__,__LINE__);
            return EXIT_SUCCESS;
         }

         int calculate_exchange_fields()
         {

            check_cuda_errors(__FILE__,__LINE__);

            if( !exchange_initialised) initialise_exchange();

            if( !empty_exchange) {
                //cudaMemcpy(d_spin3n, 				cu::atoms::d_x_spin, ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyDeviceToDevice);
                //cudaMemcpy(d_spin3n + ::atoms::num_atoms, 	cu::atoms::d_y_spin, ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyDeviceToDevice);
                //cudaMemcpy(d_spin3n + 2 * ::atoms::num_atoms, 	cu::atoms::d_z_spin, ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyDeviceToDevice);
                check_cuda_errors(__FILE__,__LINE__);

                // cusparseSpMV using CSR algorithm 1
                status = cusparseSpMV(  handle,
                                        CUSPARSE_OPERATION_NON_TRANSPOSE,
                                        &alpha,
                                        J_matrix_d,
                                        vecX,
                                        &beta,
                                        vecY,
                                        CUSPARSE_REAL,
                                        CUSPARSE_CSRMV_ALG1,
                                        spmv_buffer_d);

                if (status != CUSPARSE_STATUS_SUCCESS) {
                    std::cerr << "Matrix-vector multiplication failed" << std::endl;
                    return 1;
                }

                check_cuda_errors(__FILE__,__LINE__);

                //cudaMemcpy(cu::d_x_spin_field, d_field3n, 				::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyDeviceToDevice);
                //cudaMemcpy(cu::d_y_spin_field, d_field3n + ::atoms::num_atoms, 		::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyDeviceToDevice);
                //cudaMemcpy(cu::d_z_spin_field, d_field3n + 2 * ::atoms::num_atoms, 	::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyDeviceToDevice);

            }


            check_cuda_errors(__FILE__,__LINE__);
            return EXIT_SUCCESS;
         }


         // single_spin field device functions
         __device__ cu_real_t exchange_field_component(int *csr_cols, int* csr_rows, cu_real_t *vals, cu_real_t *s, const int i)
         {
             cu_real_t hx = 0.0;
             for (int n = csr_rows[i]; n < csr_rows[i+1]; n++){
                 int j = csr_cols[n];
                 cu_real_t val = vals[n];
                 hx += val*s[j];
             }

             return hx;
         }


      } // end namespace exchange

   } // end namespace internal

} // end namespace vcuda
