//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2016. All rights reserved.
//
//   Email: sj681@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"

namespace micromagnetic{


   //------------------------------------------------------------------------
   // Externally visiable variables
   //------------------------------------------------------------------------


   //boolean to determine whether the simulation is micromagnetic
   bool discretisation_micromagnetic = false;
   //boolean to determine whether the simulation wants stochastic fields
   bool stochastic = true;

   std::vector < std::vector <int > > P;
   std::vector < int > P1D;
   double mean_M;
   int counter;

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside micromagnetic module
      //------------------------------------------------------------------------

      //holds the cell parameters
      std::vector<double> A;
      std::vector<double> alpha;
      std::vector<double> chi_perp;
      std::vector<double> chi_para;
      std::vector<double> gamma;
      std::vector<double> ku;
      std::vector<double> ms;
      std::vector<double> Tc;

      //holds the normalised magnetisation in x,y,z
      std::vector<double> x_array;
      std::vector<double> y_array;
      std::vector<double> z_array;

      //external field vector
      std::vector<double> ext_field;


      //euler and heun arrays
      std::vector<double> x_euler_array;
      std::vector<double> y_euler_array;
      std::vector<double> z_euler_array;
      std::vector<double> x_heun_array;
      std::vector<double> y_heun_array;
      std::vector<double> z_heun_array;

      //where the magnetisation is stored between euler and heun array steps
      std::vector<double> mx_store;
      std::vector<double> my_store;
      std::vector<double> mz_store;

      //initial magnetisation vectors per step
      std::vector<double> mx_init;
      std::vector<double> my_init;
      std::vector<double> mz_init;

      //macrocell neighbourlists
      std::vector<double> macro_neighbour_list_start_index;
      std::vector<double> macro_neighbour_list_end_index;
      std::vector<double> macro_neighbour_list_array;



   } // end of internal namespace

} // end of micromagnetic namespace
