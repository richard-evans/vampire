//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins 2017. All rights reserved.
//
//   Email: sj681@york.ac.uk
//
//------------------------------------------------------------------------------
//
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"
#include "vutil.hpp"
#include "grains.hpp"
using namespace std;

namespace program{

   //---------------------------------------------------------------------------------------------------------------------
   // Setting algorithm for exchange bias in IrMn3
   //---------------------------------------------------------------------------------------------------------------------

   void setting_process(){
      //Error Checking.
      if(err::check==true) std::cout << "program::setting has been called" << std::endl;

      //Arrays to determine the 8 possible ground state spin orientations configurations
      double Configurations[8][3] = {
         {0,     0,      -1},
         {-0.94, 0,      0.3},
         {0.44,  -0.845, 0.3},
         {0.44,  0.845,  0.3},
         {0,     0,      1},
         {0.94,  0,      -0.3},
         {-0.44, 0.845,  -0.3},
         {-0.44, -0.845, -0.3}
      };
      // List the 8 possible combinations of spin orientation
      double Possible_Arrays[8][4] = {
         {0,1,2,3},
         {1,0,3,2},
         {2,3,0,1},
         {3,2,1,0},
         {4,5,6,7},
         {5,4,7,6},
         {6,7,4,5},
         {7,6,5,4}
      };

      std::vector< std::vector <int> > No_in_Sublattice;
      // Number sublattices in each grain
      No_in_Sublattice.resize(4);
      // resize array for correct number of atoms in each grain
      for(int i = 0; i < 4; i++) No_in_Sublattice[i].resize(grains::num_grains,0.0);
      std::vector <int> Total_Sub (grains::num_grains*4,0);
      std::vector <int> Local_Sub (grains::num_grains*4,0);
      std::vector <int> Chosen_array(grains::num_grains,0);
      std::vector <int> Largest_Sublattice(grains::num_grains,0);
      std::vector <int> Max_atoms(grains::num_grains,0);


      #ifdef MPICF
      stats::num_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
      #else
      stats::num_atoms = atoms::num_atoms;
      #endif

      std::vector <bool> atom_included(stats::num_atoms,false);

  //    std::cout << stats::num_atoms <<std::endl;

      //Calculates how many atoms are in the top layer of each sublattice in each grain.
      for (int atom = 0; atom < stats::num_atoms; atom++){

        int start = atoms::neighbour_list_start_index[atom];
        int end =  atoms::neighbour_list_end_index[atom] +1;
        int mat_i = atoms::type_array[atom];
        int grain = atoms::grain_array[atom];
        atom_included[atom] =false;
        if (mat_i < 4){
    //      std::cout << "A" <<atom << '\t' << atoms::x_coord_array[atom] << '\t' << atoms::y_coord_array[atom] << '\t' << atoms::z_coord_array[atom] << '\t' << atom_included[atom] << "\t" <<  atoms::type_array[atom] <<std::endl;

          for (int neighbour = start; neighbour < end; neighbour ++){
            if (atom_included[atom] == false){
            // explain what if statement is testing

            int mat_j = atoms::type_array[atoms::neighbour_list_array[neighbour]];
            std::cout << mat_j <<std::endl;
            if (mat_j == 4 || mat_j == 5){
        //      std:: cout << atom << '\t' << mat_j << std::endl;
              No_in_Sublattice[mat_i][grain]++;
              atom_included[atom] = true;
            }
            }
         }
       }
      }
    //  for (int atom = 0; atom < stats::num_atoms; atom++){
  //      std::cout << atom << '\t' << atoms::x_coord_array[atom] << '\t' << atoms::y_coord_array[atom] << '\t' << atoms::z_coord_array[atom] << '\t' << atom_included[atom] << "\t" <<  atoms::type_array[atom] <<std::endl;
//      }

      int k = 0;
      for (int j = 0; j < grains::num_grains; j ++){
         for (int i = 0; i < 4; i ++){
            Total_Sub[k] = No_in_Sublattice[i][j];


            k++;

         }
      }

      #ifdef MPICF
         MPI_Allreduce(MPI_IN_PLACE, &Total_Sub[0],grains::num_grains*4, MPI_INT, MPI_SUM,MPI_COMM_WORLD);
      #endif




      int l =0;


      //calculates which sublattice contains the most atoms for each grain.
      for (int j = 0; j < grains::num_grains; j++){

         for( int i = 0; i <4; i ++){
            if ((Total_Sub[l] > Max_atoms[j]) & (Total_Sub[l] != 0)){
               Largest_Sublattice[j] = i;
               Max_atoms[j] =Total_Sub[l];
            }
            l++;
         }
         //if (Total_Sub[l-1] !=0 )
         cerr << Largest_Sublattice[j] <<endl;
      }
      int j = 0;


      for (int i = 0; i < grains::num_grains*4; i= i +4){
        // if (Total_Sub[i] != 0){
            cout <<"number in each sublattice for grain" << j << ";" <<  Total_Sub[i] << "\t" << Total_Sub[i +1] << "\t" << Total_Sub[i+2] << "\t" << Total_Sub[i+3] << "\t" << std::endl;
            zlog << zTs() <<"number in each sublattice for grain" << j << ";" <<  Total_Sub[i] << "\t" << Total_Sub[i +1] << "\t" << Total_Sub[i+2] << "\t" << Total_Sub[i+3] << "\t" << std::endl;
            j++;
      //   }
      }


      double result,angle;
      double min_angle = 1000;
      double Direction_Closest_to_Field;

      //Loop over all possible combinations to calculate which possible vector is closest to the applied field.
      for (int i = 0; i < 8; i ++){
         result = Configurations[i][0]*sim::H_vec[0] + Configurations[i][1]*sim::H_vec[1] + Configurations[i][2]*sim::H_vec[2];
         angle = acos(result);
         //Calculates the minimum angle between the applied field and the spin directions.
         if ( angle< min_angle){
            Direction_Closest_to_Field = i;
            min_angle = angle;
         }
      }

      //Sets the sublattice with the largest number of atoms along the direction nearest the field
      //This minimises S.H
      for (int j = 0; j < grains::num_grains; j ++){
         for (int i = 0; i <8; i ++){
            if (Possible_Arrays[i][Largest_Sublattice[j]] == Direction_Closest_to_Field){
               Chosen_array[j] = i;
               break;
            }
         }
      }

      int Array;
      for (int i = 0; i < stats::num_atoms; i++){
         if(atoms::type_array[i] > 3){
            atoms::x_spin_array[i] = sim::H_vec[0];
            atoms::y_spin_array[i] = sim::H_vec[1];
            atoms::z_spin_array[i] = sim::H_vec[2];

         }

         else {

            int Array = Possible_Arrays[Chosen_array[atoms::grain_array[i]]][atoms::type_array[i]];
            atoms::x_spin_array[i] = Configurations[Array][0];
            atoms::y_spin_array[i] = Configurations[Array][1];
            atoms::z_spin_array[i] = Configurations[Array][2];
         }
      }

      // Calculate magnetisation statistics
      stats::mag_m();

      // Output data
      vout::data();

   }
}
