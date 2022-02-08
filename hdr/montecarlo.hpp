//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Adam Laverack and Richard Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef MONTECARLO_H_
#define MONTECARLO_H_

// C++ standard library headers
#include <string>
#include <vector>

// Vampire headers
#include "montecarlo.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions for montecarlo module
//--------------------------------------------------------------------------------
namespace montecarlo{

   //-----------------------------------------------------------------------------
   // Function to initialise montecarlo module
   //-----------------------------------------------------------------------------
   void initialize();

   //---------------------------------------------------------------------------
   // Function to process input file parameters for montecarlo module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

   //---------------------------------------------------------------------------
   // Function to perform one monte carlo, constrained monte carlo, or hybrid
   // cmc-mc step, respectively
   //---------------------------------------------------------------------------
   void mc_step(std::vector<double> &x_spin_array, std::vector<double> &y_spin_array, std::vector<double> &z_spin_array, int num_atoms, std::vector<int> &type_array);
   int cmc_step();
   int cmc_mc_step();
   void mc_step_parallel(std::vector<double> &x_spin_array, std::vector<double> &y_spin_array, std::vector<double> &z_spin_array, std::vector<int> &type_array);

   //---------------------------------------------------------------------------
   // Provide access to CMCinit and CMCMCinit for cmc_anisotropy and
   // hybrid_cmc programs respectively
   //---------------------------------------------------------------------------
   void CMCinit();
   void CMCMCinit();
   void mc_parallel_init(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z,
                         double min_dim[3], double max_dim[3]);

   extern bool mc_parallel_initialized;

   //---------------------------------------------------------------------------
   // Function to perform monte carlo preconditioning
   //---------------------------------------------------------------------------
   void monte_carlo_preconditioning();

   //---------------------------------------------------------------------------
   // CMC namespace; must be given external access to be initialised by
   // match.cpp and, and data can be accessed by sim.cpp and outputfunctions.cpp
   //---------------------------------------------------------------------------

   namespace cmc{

   	class cmc_material_t {
   	public:

   		double constraint_phi; /// Constrained minimisation vector (azimuthal) [degrees]
   		double constraint_phi_min; /// loop angle min [degrees]
   		double constraint_phi_max; /// loop angle max [degrees]
   		double constraint_phi_delta; /// loop angle delta [degrees]

   		double constraint_theta; /// Constrained minimisation vector (rotational) [degrees]
   		double constraint_theta_min; /// loop angle min [degrees]
   		double constraint_theta_max; /// loop angle max [degrees]
   		double constraint_theta_delta; /// loop angle delta [degrees]

   		// performance optimised rotational matrices
   		double ppolar_vector[3];
   		double ppolar_matrix[3][3];
   		double ppolar_matrix_tp[3][3];

   		// vector magnetisation
   		double M_other[3];

      	cmc_material_t():
      		constraint_phi(0.0),
      		constraint_phi_min(0.0),
      		constraint_phi_max(0.0),
      		constraint_phi_delta(5.0),
      		constraint_theta(0.0),
      		constraint_theta_min(0.0),
      		constraint_theta_max(0.0),
      		constraint_theta_delta(5.0)

      	{

      	//for(int i=0;i<100;i++){
      	//	geometry_coords[i][0]=0.0;
      	//	geometry_coords[i][1]=0.0;
      	//}
      }

   	};

   	extern std::vector<cmc_material_t> cmc_mat;

   	extern bool is_initialised;

   	extern int active_material; /// material in current hybrid loop

   	extern std::vector<std::vector< int > > atom_list;
   	extern double mc_success;
   	extern double mc_total;
   	extern double sphere_reject;
   	extern double energy_reject;

   } //end of cmc namespace

} // end of montecarlo namespace

#endif //MONTECARLO_H_
