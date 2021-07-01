//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins, Andrea Meo and Richard F L Evans 2017.
//       All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>


// Vampire headers
#include "dipole.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "cells.hpp"
#include "errors.hpp"
#include "sim.hpp"
#include "atoms.hpp"
#include "vutil.hpp"


#ifdef FFT
#include <fftw3.h>
#endif

#ifdef FFTW_OMP
#include <omp.h>
#endif

// dipole module headers
#include "internal.hpp"

// shorthand form for dipole module
namespace dp = dipole::internal;

namespace dipole{

    namespace internal{


        //------------------------------------------------------------------------------
        // Calculation of the dipole field using a Fast Fourier Transform
        //------------------------------------------------------------------------------
        // The demagnetizing field of a body can be expressed as
        //
        //                                   H = - N . M
        //
        // where H is the dipole field, N is the demagnetizing tensor, and M is the
        // magnetization. Applying discrete fourier transforms (DFT) gives
        //
        //                            DFT(H) = -DFT(N) . DFT(M)
        //
        // By rearrangement we recover the demagnetizing field H:
        //
        //                          H = iDFT[ -DFT(N) . DFT(M) ]
        //
        // Where iDFT is an inverse discrtete fourier transform. If the interaction
        // tensor N can be assumed to be translationally invariant then we can
        // accelerate the calculation of the field with a Fast Fourier Transform (FFT):
        //
        //                          H = iFFT[ -FFT(N) . FFT(M) ]
        //
        // where FFT(N) is fixed and scales with n log n numercial operations where n
        // is the number of cells.
        //
        //------------------------------------------------------------------------------

#ifdef FFT
        bool FFT_initialised = false;

        double          *M_r;       // Spatial Magnetisation
        double          *H_r;    // Spatial dipole field
        fftw_complex    *int_mat_k; // K-space interaction matrix
        fftw_complex    *M_k;       // K-space Magnetisation
        fftw_complex    *H_k;       // K-space field

        fftw_plan       plan_M;
        fftw_plan       plan_H;

        int N;
        int Nx, Ny, Nz;
        int Ncells_x, Ncells_y, Ncells_z;

        double dx, dy, dz;

        int count = 0;
        double avg_time = 0;

        std::vector<int>    cell_idx;

        inline double PBC ( double dx, double L, bool bounds) { return (bounds) ? dx - floor( (dx/L) + 0.5) * L : dx;}
        inline void complex_multiply_add( fftw_complex& a, fftw_complex& b, fftw_complex& c)
        {
            a[0] += (b[0] * c[0] - b[1] * c[1]);
            a[1] += (b[0] * c[1] + b[1] * c[0]);
        }

        double delta_fn ( const int i, const int j) { return (i==j) ? 1 : 0;}


        //-----------------------------------------------------------------------------
        // Function to initialise dipole field calculation using FFT solver
        //-----------------------------------------------------------------------------
#endif
        void initialize_fft_solver(){

#ifdef FFT

            if( fftw_init_threads() == 0)
                std::cout << "Error initialising threads for FFTW!" << std::endl;

            int Nthreads = 1;

#ifdef FFTW_OMP
            Nthreads = omp_get_max_threads();
            std::cout << "Planning FFT with Nthreads = " << Nthreads << std::endl;
#endif
            fftw_plan_with_nthreads(Nthreads);

            const double prefactor = 0.9274009994; // mu_0 * muB / (4*pi*Angstrom^3) = 1.0e-7 * 9.274009994e-24 / 1.0e-30 = 0.9274009994

            double Lx = cs::system_dimensions[0];
            double Ly = cs::system_dimensions[1];
            double Lz = cs::system_dimensions[2];


            // determine number of cells in x and y (global)
            Ncells_x = static_cast<unsigned int>(ceil((cs::system_dimensions[0]+0.01)/cells::macro_cell_size_x));
            Ncells_y = static_cast<unsigned int>(ceil((cs::system_dimensions[1]+0.01)/cells::macro_cell_size_y));
            Ncells_z = static_cast<unsigned int>(ceil((cs::system_dimensions[2]+0.01)/cells::macro_cell_size_z));

            Nx = Ncells_x;
            Ny = Ncells_y;
            Nz = Ncells_z;


            // Calculate the discretisation of each mesh point
            dx = Lx/Nx;
            dy = Ly/Ny;
            dz = Lz/Nz;

            //std::cerr << "System size is " << Lx << "  " << Ly << " " << Lz << std::endl;
            //std::cerr << "Discretisation is " << dx << "  " << dy << " " << dz << std::endl;

            // If the system is not periodic in each direction then pad
            if( !cs::pbc[0] ) Nx *=2;
            if( !cs::pbc[1] ) Ny *=2;
            if( !cs::pbc[2] ) Nz *=2;


            N = Nx*Ny*Nz;

            // Allocate 4d arrays for the magnetisation and field
            // Real space
            M_r = (double*) fftw_malloc( sizeof(double) * 3 * N);
            H_r = (double*) fftw_malloc( sizeof(double) * 3 * N);

            // complex K-space
            M_k = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * 3 * N);
            H_k = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * 3 * N);


            // The interaction matrix 4d array in K-space
            int_mat_k = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * 9 * N);

            // Calculate memory requirements and inform user
            const double mem = double(N) * ( 6*sizeof(double)  + 15*sizeof(fftw_complex)) / 1.0e6;
            zlog << zTs() << "Macrocell FFT dipole field calculation has been enabled and requires " << mem << " MB of RAM" << std::endl;
            std::cout     << "Macrocell FFT dipole field calculation has been enabled and requires " << mem << " MB of RAM" << std::endl;


            // create FFTW plans to act on the M and H arrays
            int n[] = {Nx, Ny, Nz};
            int rank = 3;
            int howmany = 3;
            int idist = 1, odist = 1;
            int istride = 3, ostride = 3;
            int *inembed = n, *onembed = n;

            // Here we plan the transforms, making use of the inter-leaved memory
            // From real space to K-space uses a real to complex
            plan_M = fftw_plan_many_dft_r2c( rank, n, howmany,
                    M_r, inembed, istride, idist,
                    M_k, onembed, ostride, odist,
                    FFTW_MEASURE);

            // From K-space to real space uses a complex to real
            plan_H = fftw_plan_many_dft_c2r( rank, n, howmany,
                    H_k, inembed, istride, idist,
                    H_r, onembed, ostride, odist,
                    FFTW_MEASURE);


            // Now setup the interaction matrix
            // Create real space array
            double *int_mat_r;
            int_mat_r = (double*) fftw_malloc( sizeof(double) * 9 * N);

            howmany = 9;
            istride = 9;
            ostride = 9;
            fftw_plan p;

            p = fftw_plan_many_dft_r2c( rank, n, howmany,
                    int_mat_r, inembed, istride, idist,
                    int_mat_k, onembed, ostride, odist,
                    FFTW_ESTIMATE);

            // Fill the interaction matrix array
            // loop over the system mesh to
            // construct the interaction matrix
            // w(r) = (\mu_0 / 4 pi r^5) (3 r \outer r - I r*r)
            int ind = 0;
            for( int i = 0 ; i < Nx; i++) {
                for( int j = 0; j < Ny; j++) {
                    for( int k = 0; k < Nz; k++) {
                            double rx = ( i > Nx/2) ? i - Nx : i;
                            double ry = ( j > Ny/2) ? j - Ny : j;
                            double rz = ( k > Nz/2) ? k - Nz : k;
                            double rij[3] = {rx*dx, ry*dy, rz*dz};
                            double r = sqrt( rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
                            double w0 = prefactor / (r*r*r*r*r);
                            // FFTW does not normalise the transform so we do here.
                            w0 /= double(N);
                            for( int a = 0; a < 3; a++) {
                                for( int b = 0; b < 3; b++) {
                                    int id = b + 3*a;
                                    int idx = id + 9 * (k + Nz * (j + Ny * i));
                                    int_mat_r[idx] = w0 * ( 3* rij[a] * rij[b] - delta_fn(a,b) * r * r);
                                }
                            }
                        }
                }
            }
            // Zero out the self-interaction
            for( int a = 0; a < 9; a++)
                int_mat_r[a] = 0.0;

            /*
            for( int i = 0 ; i < Nx; i++) {
                for( int j = 0; j < Ny; j++) {
                    for( int k = 0; k < Nz; k++) {
                        std::cerr << i << "  " << j << "  " << k << "  ";
                            for( int a = 0; a < 3; a++) {
                                for( int b = 0; b < 3; b++) {
                                    int id = b + 3*a;
                                    int idx = id + 9 * (k + Nz * (j + Ny * i));
                                    std::cerr << int_mat_r[idx] << "  ";
                                }
                            }
                            std::cerr << std::endl;
                        }
                }
            } */

            // Now perform the FFT to get the K-space values
            fftw_execute(p);

            // Plan is no longer need so can be destroyed
            fftw_destroy_plan(p);

            fftw_free(int_mat_r);


            // Now construct list of where each cell sits in the arrays
            cell_idx.resize(cells::num_cells);

            for ( int i=0; i < Ncells_x; i++){
                for ( int j=0; j < Ncells_y; j++){
                    for ( int k=0; k < Ncells_z; k++){
                        //std::cerr << i << "  " << rx << "  " << ry << "  " << rz << "  " << ix << "  " << iy << "  " << iz << std::endl;
                        cell_idx[i] = k + Nz * (j + Ny * i);
                    }
                }
            }

            //Allocate storage for the cell field
            dipole::cells_field_array_x.resize(cells_num_cells,0.0);
            dipole::cells_field_array_y.resize(cells_num_cells,0.0);
            dipole::cells_field_array_z.resize(cells_num_cells,0.0);

            // resize mu_0*Hd-field cells array
            dipole::cells_mu0Hd_field_array_x.resize(cells_num_cells,0.0);
            dipole::cells_mu0Hd_field_array_y.resize(cells_num_cells,0.0);
            dipole::cells_mu0Hd_field_array_z.resize(cells_num_cells,0.0);

            FFT_initialised = true;

#endif

        }

        void update_field_fft(){


#ifdef FFT
            if(!FFT_initialised) {
                std::cout << "FFT dipole has been called but not initialised." << std::endl;
                exit(-1);
            }

            // instantiate timer
            vutil::vtimer_t timer;


            // update cell magnetisations
            cells::mag();

            //   start timer
            timer.start();

            for ( int i = 0; i < 3*N; i++)
                M_r[i] = 0.0;

            for( int cell = 0; cell < cells::num_cells; cell++){
                int c_idx = cell_idx[cell];
                const double mu = 1.0; //atoms::m_spin_array[atom];
                M_r[3*c_idx]     = cells::mag_array_x[cell]*mu;
                M_r[3*c_idx + 1] = cells::mag_array_y[cell]*mu;
                M_r[3*c_idx + 2] = cells::mag_array_z[cell]*mu;
            }



            // Forward FFT to get M_k
            fftw_execute(plan_M);


            // H_k is the product of int_mat_k and M_k
            for( int i = 0 ; i < Nx; i++) {
                for( int j = 0; j < Ny; j++) {
                    for( int k = 0; k < Nz; k++) {
                        for( int a = 0; a < 3; a++) {
                            for( int b = 0; b < 3; b++) {
                                int id = b + 3*a;
                                int int_idx = id + 9 * (k + Nz * (j + Ny * i));
                                int M_idx = a + 3*(k + Nz * (j + Ny * i));
                                int H_idx = b + 3*(k + Nz * (j + Ny * i));
                                complex_multiply_add( H_k[H_idx], int_mat_k[int_idx], M_k[M_idx] );
                            }
                        }
                    }
                }
            }

            // Inverse FFT to get H in real space
            fftw_execute( plan_H );

            // save total dipole field to atomic field array
            for ( int cell = 0; cell < cells::num_cells; cell++) {
                int c_idx = cell_idx[cell];
                dipole::cells_field_array_x[cell] = H_r[3*c_idx];
                dipole::cells_field_array_y[cell] = H_r[3*c_idx + 1];
                dipole::cells_field_array_z[cell] = H_r[3*c_idx + 2];
            }
            // For MPI version, only add local atoms
#ifdef MPICF
            const int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
#else
            const int num_local_atoms = dipole::internal::num_atoms;
#endif

            // Update Atomistic Dipolar Field and Demag Field Array
            for(int atom=0;atom<num_local_atoms;atom++){

                const int cell = dipole::internal::atom_cell_id_array[atom];

                int type = dipole::internal::atom_type_array[atom];

                if(dipole::internal::cells_num_atoms_in_cell[cell]>0 && mp::material[type].non_magnetic==false){

                    // Copy B-field from macrocell to atomistic spin
                    // Copy B-field from macrocell to atomistic spin
                    dipole::atom_dipolar_field_array_x[atom] = dipole::cells_field_array_x[cell];
                    dipole::atom_dipolar_field_array_y[atom] = dipole::cells_field_array_y[cell];
                    dipole::atom_dipolar_field_array_z[atom] = dipole::cells_field_array_z[cell];

                    // Unroll Hdemag field
                    //dipole::atom_mu0demag_field_array_x[atom] = dipole::cells_mu0Hd_field_array_x[cell];
                    //dipole::atom_mu0demag_field_array_y[atom] = dipole::cells_mu0Hd_field_array_y[cell];
                    //dipole::atom_mu0demag_field_array_z[atom] = dipole::cells_mu0Hd_field_array_z[cell];

                }
            }

            timer.stop();

            avg_time += timer.elapsed_time();
            count++;

#endif

            return;

        }

        //-----------------------------------------------------------------------------
        // Function to finalize FFT solver and release memory
        //-----------------------------------------------------------------------------
        void finalize_fft_solver(){
#ifdef FFT
            std::cout << "Average FFT dipole compute time = " << avg_time / double(count) << std::endl;
            zlog << zTs() << "Average FFT dipole compute time = " << avg_time / double(count) << std::endl;

            // Print informative message to log file and screen
            std::cout << "Deallocating memory for FFT dipole calculation" << std::endl;
            zlog << zTs() << "Deallocating memory for FFT dipole calculation" << std::endl;

            // Free memory from FFT complex variables
            fftw_free(M_r);
            fftw_free(M_k);
            fftw_free(H_r);
            fftw_free(H_k);
            fftw_free(int_mat_k);

            fftw_destroy_plan(plan_M);
            fftw_destroy_plan(plan_H);

            fftw_cleanup_threads();
#endif
            return;

        }
    }
}
