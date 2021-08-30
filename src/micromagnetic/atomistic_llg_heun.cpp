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
#include <cmath>
#include <cstdlib>
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "LLG.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "micromagnetic.hpp"
#include "internal.hpp"

namespace mm = micromagnetic::internal;


namespace atomistic_LLG_arrays{

   // Local arrays for LLG integration
   std::vector <double> x_euler_array;
   std::vector <double> y_euler_array;
   std::vector <double> z_euler_array;

   std::vector <double> x_heun_array;
   std::vector <double> y_heun_array;
   std::vector <double> z_heun_array;

   std::vector <double> x_spin_storage_array;
   std::vector <double> y_spin_storage_array;
   std::vector <double> z_spin_storage_array;

   std::vector <double> x_initial_spin_array;
   std::vector <double> y_initial_spin_array;
   std::vector <double> z_initial_spin_array;

   bool LLG_set=false; ///< Flag to define state of LLG arrays (initialised/uninitialised)

}


namespace micromagnetic{


   int atomistic_LLGinit(){

      // check calling of routine if error checking is activated
      if(err::check==true){std::cout << "micromagnetic::atomistic_LLG_init has been called" << std::endl;}

      using namespace atomistic_LLG_arrays;
      //resize LLG arrays
      x_spin_storage_array.resize(atoms::num_atoms,0.0);
      y_spin_storage_array.resize(atoms::num_atoms,0.0);
      z_spin_storage_array.resize(atoms::num_atoms,0.0);

      x_initial_spin_array.resize(atoms::num_atoms,0.0);
      y_initial_spin_array.resize(atoms::num_atoms,0.0);
      z_initial_spin_array.resize(atoms::num_atoms,0.0);

      x_euler_array.resize(atoms::num_atoms,0.0);
      y_euler_array.resize(atoms::num_atoms,0.0);
      z_euler_array.resize(atoms::num_atoms,0.0);

      x_heun_array.resize(atoms::num_atoms,0.0);
      y_heun_array.resize(atoms::num_atoms,0.0);
      z_heun_array.resize(atoms::num_atoms,0.0);

      LLG_set=true;

      return EXIT_SUCCESS;
   }


   int atomistic_LLG_Heun(){


      // check calling of routine if error checking is activated
      if(err::check==true){std::cout << "micromagnetic::LLG_Heun has been called" << std::endl;}

      using namespace atomistic_LLG_arrays;

      // Check for initialisation of LLG integration arrays
      if(LLG_set== false) atomistic_LLGinit();

      // Local variables for system integration
      //const int num_atoms=atoms::num_atoms;
      double xyz[3];		// Local Delta Spin Components
      double S_new[3];	// New Local Spin Moment
      double mod_S;		// magnitude of spin moment

      // Store initial spin positions
      for(int atom_list=0;atom_list<number_of_atomistic_atoms;atom_list++){
         int atom = list_of_atomistic_atoms[atom_list];
         x_initial_spin_array[atom] = atoms::x_spin_array[atom];
         y_initial_spin_array[atom] = atoms::y_spin_array[atom];
         z_initial_spin_array[atom] = atoms::z_spin_array[atom];
      }

      //loops over sections of atoms which are numerically adjacent ie. 4,5,6,7,8,9 then 12,13,14,15
      //num calcualtions is the number of sections of adjacent atoms.
      int num_calculations = mm::fields_neighbouring_atoms_begin.size();

      //caclualtes the spin fields and external fields for each section.
      for (int i = 0; i < num_calculations; i ++){
         int begin = mm::fields_neighbouring_atoms_begin[i];		//the start and end atom index for each section.
         int end = mm::fields_neighbouring_atoms_end[i];
         sim::calculate_spin_fields(begin,end);								//calls the function from the usual atomsitic vampire
         sim::calculate_external_fields(begin,end);						//external fields on each atom from the atomistic vampire.
      }

//std::cin.get();

      //Calculate Euler gradients
      for(int atom_list=0;atom_list<number_of_atomistic_atoms;atom_list++){
         //calcualtes the atom if from the atom list
         int atom = list_of_atomistic_atoms[atom_list];
         //sets the material for each atom
         const int imaterial=atoms::type_array[atom];

         //alpha is used as 1/(1+a^2) and a/(1+a^2) so these are stored as variables.
         const double one_oneplusalpha_sq = mp::material[imaterial].one_oneplusalpha_sq; // material specific alpha and gamma
         const double alpha_oneplusalpha_sq = mp::material[imaterial].alpha_oneplusalpha_sq;

         // Store local spin in S and local field in H
         const double S[3] = {atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};
         const double H[3] = {atoms::x_total_spin_field_array[atom]+atoms::x_total_external_field_array[atom],
            atoms::y_total_spin_field_array[atom]+atoms::y_total_external_field_array[atom],
            atoms::z_total_spin_field_array[atom]+atoms::z_total_external_field_array[atom]};
            // Calculate Delta S
            xyz[0]=(one_oneplusalpha_sq)*(S[1]*H[2]-S[2]*H[1]) + (alpha_oneplusalpha_sq)*(S[1]*(S[0]*H[1]-S[1]*H[0])-S[2]*(S[2]*H[0]-S[0]*H[2]));
            xyz[1]=(one_oneplusalpha_sq)*(S[2]*H[0]-S[0]*H[2]) + (alpha_oneplusalpha_sq)*(S[2]*(S[1]*H[2]-S[2]*H[1])-S[0]*(S[0]*H[1]-S[1]*H[0]));
            xyz[2]=(one_oneplusalpha_sq)*(S[0]*H[1]-S[1]*H[0]) + (alpha_oneplusalpha_sq)*(S[0]*(S[2]*H[0]-S[0]*H[2])-S[1]*(S[1]*H[2]-S[2]*H[1]));

            // Store dS in euler array
            x_euler_array[atom]=xyz[0];
            y_euler_array[atom]=xyz[1];
            z_euler_array[atom]=xyz[2];

            // Calculate Euler Step
            S_new[0]=S[0]+xyz[0]*mp::dt;
            S_new[1]=S[1]+xyz[1]*mp::dt;
            S_new[2]=S[2]+xyz[2]*mp::dt;

            // Normalise Spin Length
            mod_S = 1.0/sqrt(S_new[0]*S_new[0] + S_new[1]*S_new[1] + S_new[2]*S_new[2]);
            S_new[0]=S_new[0]*mod_S;
            S_new[1]=S_new[1]*mod_S;
            S_new[2]=S_new[2]*mod_S;

            //Writing of Spin Values to Storage Array
            x_spin_storage_array[atom]=S_new[0];
            y_spin_storage_array[atom]=S_new[1];
            z_spin_storage_array[atom]=S_new[2];
         }

         // Copy new spins to spin array
         //Calculate Euler gradients
         for(int atom_list=0;atom_list<number_of_atomistic_atoms;atom_list++){
            //calcualtes the atom if from the atom list
            int atom = list_of_atomistic_atoms[atom_list];
       		atoms::x_spin_array[atom]=x_spin_storage_array[atom];
       		atoms::y_spin_array[atom]=y_spin_storage_array[atom];
       		atoms::z_spin_array[atom]=z_spin_storage_array[atom];
       	}



         //recalcualtes the spin fields for each section as before
         for (int i = 0; i < num_calculations; i ++){
            int begin = mm::fields_neighbouring_atoms_begin[i];
            int end = mm::fields_neighbouring_atoms_end[i];
            sim::calculate_spin_fields(begin,end);
         }


         // Calculate Heun Gradients
         for(int atom_list=0;atom_list<number_of_atomistic_atoms;atom_list++){
            int atom = list_of_atomistic_atoms[atom_list];

            //calcualtes the material and 1/(1+a^2) and a/(1+a^2)
            const int imaterial=atoms::type_array[atom];;
            const double one_oneplusalpha_sq = mp::material[imaterial].one_oneplusalpha_sq;
            const double alpha_oneplusalpha_sq = mp::material[imaterial].alpha_oneplusalpha_sq;



            // Store local spin in Sand local field in H
            const double S[3] = {atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};
            const double H[3] = {atoms::x_total_spin_field_array[atom]+atoms::x_total_external_field_array[atom] + mp::material[imaterial].pinning_field_unit_vector[0],
               atoms::y_total_spin_field_array[atom]+atoms::y_total_external_field_array[atom] + mp::material[imaterial].pinning_field_unit_vector[1],
               atoms::z_total_spin_field_array[atom]+atoms::z_total_external_field_array[atom] + mp::material[imaterial].pinning_field_unit_vector[2]};

               // Calculate Delta S
               xyz[0]=(one_oneplusalpha_sq)*(S[1]*H[2]-S[2]*H[1]) + (alpha_oneplusalpha_sq)*(S[1]*(S[0]*H[1]-S[1]*H[0])-S[2]*(S[2]*H[0]-S[0]*H[2]));
               xyz[1]=(one_oneplusalpha_sq)*(S[2]*H[0]-S[0]*H[2]) + (alpha_oneplusalpha_sq)*(S[2]*(S[1]*H[2]-S[2]*H[1])-S[0]*(S[0]*H[1]-S[1]*H[0]));
               xyz[2]=(one_oneplusalpha_sq)*(S[0]*H[1]-S[1]*H[0]) + (alpha_oneplusalpha_sq)*(S[0]*(S[2]*H[0]-S[0]*H[2])-S[1]*(S[1]*H[2]-S[2]*H[1]));

               // Store dS in heun array
               x_heun_array[atom]=xyz[0];
               y_heun_array[atom]=xyz[1];
               z_heun_array[atom]=xyz[2];
            }


            //calcualtes the new spin arrays from the euler and heun Gradients
            //S = S_initial + 1/2dt ( euler + heun)
            for(int atom_list=0;atom_list<number_of_atomistic_atoms;atom_list++){
               int atom = list_of_atomistic_atoms[atom_list];
               S_new[0]=x_initial_spin_array[atom]+mp::half_dt*(x_euler_array[atom]+x_heun_array[atom]);
               S_new[1]=y_initial_spin_array[atom]+mp::half_dt*(y_euler_array[atom]+y_heun_array[atom]);
               S_new[2]=z_initial_spin_array[atom]+mp::half_dt*(z_euler_array[atom]+z_heun_array[atom]);

               // Normalise Spin Length
               mod_S = 1.0/sqrt(S_new[0]*S_new[0] + S_new[1]*S_new[1] + S_new[2]*S_new[2]);

               S_new[0]=S_new[0]*mod_S;
               S_new[1]=S_new[1]*mod_S;
               S_new[2]=S_new[2]*mod_S;

               // Copy new spins to spin array
               atoms::x_spin_array[atom]=S_new[0];
               atoms::y_spin_array[atom]=S_new[1];
               atoms::z_spin_array[atom]=S_new[2];
            }
            if (enable_resistance && mm::resistance_layer_2 != mm::resistance_layer_1)  micromagnetic::MR_resistance = mm::calculate_resistance();

            return EXIT_SUCCESS;
         }

      }
