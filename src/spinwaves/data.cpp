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

// C++ standard library headers

// Vampire headers
#include "spinwaves.hpp"

// sw module headers
#include "internal.hpp"

namespace spinwaves{

   //------------------------------------------------------------------------------
   // Externally visible variables ------------------------------------------------
   //------------------------------------------------------------------------------
   std::vector <double> skx_r;
   std::vector <double> skx_i;
   std::vector <double> skx_r_node;
   std::vector <double> skx_i_node;
   std::vector <double> skx_r_node_transposed;
   std::vector <double> skx_i_node_transposed;

   // THIS IS TEMPORARY
   int nk_per_rank;
   int scatterlength;
   std::vector<double> skx_r_scatter;
   std::vector<double> skx_i_scatter;

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside sw module
      //------------------------------------------------------------------------
      bool enabled; // bool to enable module
      std::vector <double> kx;
      std::vector <double> ky;
      std::vector <double> kz;
      std::vector<internal::mp_t> mp; // array of material properties

      // JRH Internally visible path vectors
      std::vector<double> pathx;
      std::vector<double> pathy;
      std::vector<double> pathz;

      // JRH internally visible filename
      std::string filename;
      std::string filetype = "path";

      // fft-in-time options;
      std::string reduc_ver = "rank0";
      std::vector<std::string> component;
      std::vector<std::vector<double>*> sw_array;


      bool isf;
      std::vector<bool> oss;
      std::vector<bool> cm;
      std::vector<bool> normk;
      std::vector<int> mat_in_spec;
      std::vector<int> mat;
      int nspec = 1;
      bool prefactor = true;


      std::vector<double> cos_k;
      std::vector<double> sin_k;
      std::vector<int> atom_mask;
      std::vector<int> spec_mask;


      // JRH number of time and kpoints
      int nk;
      int nt;

      // to check number of spectrums
      std::vector<int> super_index_values;


   } // end of internal namespace

} // end of sw namespace
