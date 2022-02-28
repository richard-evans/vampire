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
      //    std::cout << grains::num_grains <<std::endl;
      for(int i = 0; i < 4; i++) No_in_Sublattice[i].resize(grains::num_grains,0.0);

      std::vector <int> Local_Sub (grains::num_grains*4,0);
      std::vector <int> Chosen_array(grains::num_grains,0);
      std::vector <int> Largest_Sublattice(grains::num_grains,0);
      std::vector <int> Max_atoms(grains::num_grains,0);

      //#ifdef MPICF
      //   atoms::num_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
      //#else
      //   atoms::num_atoms = atoms::num_atoms;
      //#endif

      //std::cerr << "A" << atoms::num_atoms << "\t" << grains::num_grains << std::endl;

      //Calculates how many atoms are in the top layer of each sublattice in each grain.
      for (int atom = 0; atom < atoms::num_atoms; atom++){

         for (int neighbour = atoms::neighbour_list_start_index[atom]; neighbour < atoms::neighbour_list_end_index[atom]; neighbour ++){
            // explain what if statement is testing - yes Sarah...
            //std::cout << atom << "\t" << neighbour <<atoms::type_array[atom] << '\t' << atoms::type_array[atoms::neighbour_list_array[neighbour]] << std::endl;
            if ((atoms::type_array[atom] >3) && atoms::type_array[atom] != 8 && (atoms::type_array[atoms::neighbour_list_array[neighbour]] < 4) && atoms::type_array[atoms::neighbour_list_array[neighbour]] != 8){
               //         std::cerr << atoms::grain_array[atom] << "\t"<< atoms::type_array[atoms::neighbour_list_array[neighbour]] << "\t" << atoms::type_array[atom] <<endl;
               No_in_Sublattice[atoms::type_array[atoms::neighbour_list_array[neighbour]]][atoms::grain_array[atom]]++;
               // cerr << No_in_Sublattice[atoms::type_array[atoms::neighbour_list_array[neighbour]]][atoms::grain_array[atom]]<<endl;
            }
         }
      }

      int k = 0;
      for (int j = 0; j < grains::num_grains; j ++){
         for (int i = 0; i < 4; i ++){
            Local_Sub[k] = No_in_Sublattice[i][j];
            k++;

         }

      }
   //   std::cerr << "before" << Local_Sub[0] << '\t' << Local_Sub[1] << '\t' << Local_Sub[2] << '\t' << Local_Sub[3] << std::endl;

      #ifdef MPICF
         MPI_Allreduce(MPI_IN_PLACE, &Local_Sub[0],grains::num_grains*4,MPI_INT,MPI_SUM, MPI_COMM_WORLD);
      #endif

      //std::cerr<< "after" << Local_Sub[0] << '\t' << Local_Sub[1] << '\t' << Local_Sub[2] << '\t' << Local_Sub[3] << std::endl;

      int l =0;

      //calculates which sublattice contains the most atoms for each grain.
      for (int j = 0; j < grains::num_grains; j++){

         for( int i = 0; i <4; i ++){
            //       std::cout << Local_Sub[l] << "\t" << Local_Sub[l] << std::endl;
            if ((Local_Sub[l] > Max_atoms[j]) & (Local_Sub[l] != 0)){
               Largest_Sublattice[j] = i;
               Max_atoms[j] =Local_Sub[l];
            }
            l++;
         }
         if (Local_Sub[l-1] !=0 ) cerr << Largest_Sublattice[j] <<endl;
      }
      int j = 0;


      for (int i = 0; i < grains::num_grains*4; i= i +4){
         if (Local_Sub[i] != 0){
            cerr <<"number in each sublattice for grain" << j << ";" <<  Local_Sub[i] << "\t" << Local_Sub[i +1] << "\t" << Local_Sub[i+2] << "\t" << Local_Sub[i+3] << "\t" << std::endl;
            zlog << zTs() <<"number in each sublattice for grain" << j << ";" <<  Local_Sub[i] << "\t" << Local_Sub[i +1] << "\t" << Local_Sub[i+2] << "\t" << Local_Sub[i+3] << "\t" << std::endl;
            j++;
         }
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
      // std::cout << Direction_Closest_to_Field <<std::endl;

      //Sets the sublattice with the largest number of atoms along the direction nearest the field
      //This minimises S.H
      for (int j = 0; j < grains::num_grains; j ++){
         for (int i = 0; i <8; i ++){
            if (Possible_Arrays[i][Largest_Sublattice[j]] == Direction_Closest_to_Field){
               //   std::cout << i <<std::endl;
               Chosen_array[j] = i;
               break;
            }
         }
      }


      for (int i = 0; i <atoms::num_atoms; i++){
         //	std::cout << atoms::type_array[i] << "\t" << atoms::z_spin_array[i] <<std::endl;
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
      stats::update();

      // Output data
      vout::data();

   }
}
