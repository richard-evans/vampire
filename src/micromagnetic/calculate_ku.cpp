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
// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"
#include "material.hpp"
#include "cells.hpp"
#include "random.hpp"
#include "../anisotropy/internal.hpp"


void crossProduct(std::vector<double> vect_A, std::vector < std::vector < double > > vect_B, std::vector < std::vector < double > > &cross_P,int cell)
{

    cross_P[cell][0] = vect_A[1] * vect_B[cell][2] - vect_A[2] * vect_B[cell][1];
    cross_P[cell][1] = vect_A[0] * vect_B[cell][2] - vect_A[2] * vect_B[cell][0];
    cross_P[cell][2] = vect_A[0] * vect_B[cell][1] - vect_A[1] * vect_B[cell][0];
}

void crossProduct2(std::vector < std::vector < double > > vect_A, std::vector < std::vector < double > > vect_B, std::vector < std::vector < double > > &cross_P,int cell)
{

    cross_P[cell][0] = vect_A[cell][1] * vect_B[cell][2] - vect_A[cell][2] * vect_B[cell][1];
    cross_P[cell][1] = vect_A[cell][0] * vect_B[cell][2] - vect_A[cell][2] * vect_B[cell][0];
    cross_P[cell][2] = vect_A[cell][0] * vect_B[cell][1] - vect_A[cell][1] * vect_B[cell][0];
}

namespace micromagnetic{
   namespace internal{

      //calculates the magnetocrystaline anisotropy for each cell from the individual atomic anisotropies
      //ku = sum (ku) for each cell
      std::vector<double> calculate_ku(const int num_atoms,
                                       const int num_cells,
                                       std::vector<int> cell_array,                  //1D array storing which cell each atom is in
                                       const std::vector<int> type_array,            //1D array storing which material each atom is
                                       std::vector <mp::materials_t> material){      //class of material parameters for the atoms



         std::vector < std::vector < double > > ku_uv;
         ku_uv.resize(num_cells);
         for (int i = 0; i < num_cells; ++i)
         ku_uv[i].resize(3,0.0);

         std::vector < std::vector < double > > e1;
         e1.resize(num_cells);
         for (int i = 0; i < num_cells; ++i)
         e1[i].resize(3,0.0);

         std::vector < std::vector < double > >  e2;
         e2.resize(num_cells);
         for (int i = 0; i < num_cells; ++i)
         e2[i].resize(3,0.0);

         std::vector<double> ku(num_cells,0.0);
         std::vector<double> RN(3,0.0);


         RN[0] = mtrandom::gaussian();
         RN[1] = mtrandom::gaussian();
         RN[2] = mtrandom::gaussian();
         const double normal = sqrt(RN[0]*RN[0] + RN[1]*RN[1]  + RN[2]*RN[2]);
         RN[0] = RN[0]/normal;
         RN[1] = RN[1]/normal;
         RN[2] = RN[2]/normal;


         //sums over all atoms to add the ku per cell
         for (int atom = 0; atom < num_atoms; atom++){
            int cell = cell_array[atom];
            int mat = type_array[atom];
            ku[cell] = ku[cell] - anisotropy::internal::mp[mat].ku2;//mp::material[mat].Ku1_SI; //need to add a function here to access anisotropy module
            ku_x[cell] = ku_x[cell] + anisotropy::internal::mp[mat].ku_vector[0];
            ku_y[cell] = ku_y[cell] + anisotropy::internal::mp[mat].ku_vector[1];
            ku_z[cell] = ku_z[cell] + anisotropy::internal::mp[mat].ku_vector[2];
         }



         // #ifdef MPICF
         //    MPI_Allreduce(MPI_IN_PLACE, &ku_x[0],     num_cells,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
         //    MPI_Allreduce(MPI_IN_PLACE, &ku_y[0],     num_cells,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
         //    MPI_Allreduce(MPI_IN_PLACE, &ku_z[0],     num_cells,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
         // #endif

         for (int lc = 0; lc < cells::num_local_cells; lc++){
           int cell = cells::local_cell_array[lc];

           double normal = sqrt(ku_x[cell]*ku_x[cell] + ku_y[cell]*ku_y[cell]  + ku_z[cell]*ku_z[cell]);
           ku_uv[cell][0] =ku_x[cell]/normal;
           ku_uv[cell][1] =ku_y[cell]/normal;
           ku_uv[cell][2] =ku_z[cell]/normal;


           crossProduct(RN,ku_uv,e1,cell);
           normal = sqrt(e1[cell][0]*e1[cell][0] + e1[cell][1]*e1[cell][1]  + e1[cell][2]*e1[cell][2]);
           e1[cell][0] =e1[cell][0]/normal;
           e1[cell][1] =e1[cell][1]/normal;
           e1[cell][2] =e1[cell][2]/normal;
           crossProduct2(e1,ku_uv,e2,cell);
           normal = sqrt(e2[cell][0]*e2[cell][0] + e2[cell][1]*e2[cell][1]  + e2[cell][2]*e2[cell][2]);
           e2[cell][0] =e2[cell][0]/normal;
           e2[cell][1] =e2[cell][1]/normal;
           e2[cell][2] =e2[cell][2]/normal;

           ku_x[cell] = ku[cell] * sqrt(e1[cell][0]*e1[cell][0] + e2[cell][0]*e2[cell][0]);
           ku_y[cell] = ku[cell] * sqrt(e1[cell][1]*e1[cell][1] + e2[cell][1]*e2[cell][1]);
           ku_z[cell] = ku[cell] * sqrt(e1[cell][2]*e1[cell][2] + e2[cell][2]*e2[cell][2]);
           //std::cout << cell << '\t' <<ku[cell] << '\t' <<  ku_x[cell] << '\t' <<  ku_y[cell] << '\t' <<  ku_z[cell] << '\t' << std::endl;

        }

        #ifdef MPICF
           MPI_Allreduce(MPI_IN_PLACE, &ku_x[0],     num_cells,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
           MPI_Allreduce(MPI_IN_PLACE, &ku_y[0],     num_cells,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
           MPI_Allreduce(MPI_IN_PLACE, &ku_z[0],     num_cells,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
        #endif

        // why is this summed again?




         return ku;               //returns a
      }
   } //closes the internal namspace
}  //closes the micromagnetic namespace
