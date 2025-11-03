//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sergiu Ruta 2022. All rights reserved.
//
//   Email: sergiu.ruta@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef SW_INTERNAL_H_
#define SW_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the sw module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "spinwaves.hpp"

// sw module headers
#include "internal.hpp"
#include <fftw3.h>
#include "vector"

namespace spinwaves{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      //-----------------------------------------------------------------------------
      // internal materials class for storing material parameters
      //-----------------------------------------------------------------------------
      class mp_t{

          private:

          public:

             //------------------------------
             // material parameter variables
             //------------------------------
             double test;

             // constructor
             mp_t (const unsigned int max_materials = 100):
                test(0.0) // constructor initialisation of test variable
             {
                // constructor body for initialising more complex data/arrays
             }; // end of constructor

       }; // end of internal::mp class

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------

      extern bool enabled; // bool to enable module

      extern std::vector<internal::mp_t> mp; // array of material properties
      extern std::vector <double> kx;
      extern std::vector <double> ky;
      extern std::vector <double> kz;

      // JRH extern file variables
      extern std::string filename;
      extern std::string filetype;

      // JRH path vectors
      extern std::vector<double> pathx;
      extern std::vector<double> pathy;
      extern std::vector<double> pathz;

      // JRH fourier prefactors
      extern std::vector<double> cos_k;
      extern std::vector<double> sin_k;

      // JRH mask for calculations of spinwave for specific material type
      extern std::vector<int> atom_mask;
      extern std::vector<int> spec_mask;


      //JRH time variabls
      extern int nt;
      extern int nk;

      // JRH fft-in-time variables
      extern std::string reduc_ver;
      extern std::vector<std::string> component;
      extern std::vector<std::vector<double>*> sw_array;
      extern std::vector<bool> oss, cm, normk;
      extern std::vector<int> mat, mat_in_spec;
      extern bool prefactor;
      extern bool isf;
      extern int nspec;
      extern std::vector<int> super_index_values;


      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

      void path_sc();
      extern void path_bcc();
      extern void path_fcc();
      extern void path_hcp();
      extern void path_mn2au();
      extern void path_rocksalt();
      extern void determine_path();
      extern int gcd(int a, int b);
      extern void save_frequencies();
      extern void initialise_arrays();
      extern void calculate_material_mask();
      extern void determine_spin_component();
      extern void check_numbering_of_spectrums();
      extern void calculate_fourier_prefactor(const std::vector<double>& rx, const std::vector<double>& ry, const std::vector<double>& rz);
      extern void determine_kpoints_from_user_specific_k(const double dimx, const double dimy, const double dimz,	const double uc_x, const double uc_y, const double uc_z);
      extern void determine_kpoints_from_user_high_sym_path(const double dimx, const double dimy, const double dimz,	const double uc_x, const double uc_y, const double uc_z);


      // post analysis functions
      extern void complex_magnitude(fftw_complex *os);
      extern void one_sided_spectrum(fftw_complex *os);
      extern void write_to_file(fftw_complex *os, int k, int spec);
      extern void normalise_each_kpoint(fftw_complex *os);
      extern void write_intermediate_to_file(fftw_complex *os, int k, int spec);



      // Peak finding at some point?????/
      // extern void write_peaks_to_file();
      // extern void calculate_spectrum_peaks();


   } // end of internal namespace

} // end of sw namespace

#endif //SW_INTERNAL_H_
