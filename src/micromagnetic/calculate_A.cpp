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
#include "atoms.hpp"
#include "internal.hpp"
#include "material.hpp"
#include "create.hpp"

// micromagnetic module headers
#include <stdlib.h>
#include <vector>
#include <iostream>
#include "errors.hpp"
#include "vio.hpp"

namespace micromagnetic{

   namespace internal{

      //---------------------------------------------------------------------------------------------------
      // function to calulate the intercell exchange.
      // Each cell has an exchange interaction with every other cell in the system.
      // This is calculated as a sum over all the atoms on the surface of the cells and their interaction
      // with atoms in the other cell.

      //---------------------------------------------------------------------------------------------------
      //
      //                       A = 1/4*sum((Jij) *(x_i - x_j)^2)
      //This calculates the interaction between neighbouring cells therefore only interactions
      // where i and j are in different macrocells are summed.
      //
      // we take some components of the exchange field in here and calculate them as part of a ?
      // therefore we return
      //
      //      A = (1/2)*sum((Jij)(x_i - x_j)^2) * (V_cell/V_atomic) * (1/l_cell) * (1/Ms)
      //
      //---------------------------------------------------------------------------------------
      // Nice feature of atomic scale interactions is that all cell-cell interactions are
      // always included, even in parallel version of the code (no need to reduce across all
      // processors)
      //---------------------------------------------------------------------------------------
      // note: this solution is 2N^2 in memory usage - not scalable for large numbers of cells
      // (eg 1M cells = 64*10^12 bytes- 64 TB memory PER PROCESS(!))
      // solution is partial sum via MPI as only local interactions are needed (but gets messy).


      std::vector< double > calculate_a(int num_atoms,
                                        int num_cells,
                                        int num_local_cells,
                                        std::vector<int> cell_array,                      //1D array storing which cell each atom is in
                                        std::vector<int> neighbour_list_array,            //1D vector listing the nearest neighbours of each atom
                                        std::vector<int> neighbour_list_start_index,      //1D vector storing the start index for each atom in the neighbour_list_array
                                        std::vector<int> neighbour_list_end_index,        //1D vector storing the end index for each atom in the neighbour_list_array
                                        const std::vector<int> type_array,                //1D array storing which material each cell is
                                        std::vector <mp::materials_t> material,           //class of material parameters for the atoms
                                        std::vector <double> volume_array,                //1D array storing the volume of each cell
                                        std::vector <double> x_coord_array,
                                        std::vector <double> y_coord_array,
                                        std::vector <double> z_coord_array,
                                        double num_atoms_in_unit_cell,
                                        std::vector<int> local_cell_array){             //cell aray local to each processor



            int array_index = 0;
            std::vector< std::vector< double> > a2d; //stores the 2D exchange constants.
            a2d.resize(num_cells, std::vector<double>(num_cells,0.0));

            std::vector< std::vector< double> > num_interactions_cell; //stores the 2D exchange constants.
            num_interactions_cell.resize(num_cells, std::vector<double>(num_cells,0.0));

            // For MPI version, only add local atoms
            #ifdef MPICF
               int num_local_atoms =  vmpi::num_core_atoms+vmpi::num_bdry_atoms;
            #else
               int num_local_atoms = atoms::num_atoms;
            #endif

            std::vector<double> a;
            std::vector<double> a_test(num_cells*num_cells,0.0);
            std::vector<double> N(num_cells*num_cells,0.0);

            //calculates the atomic volume  = volume of one cell/number of atoms in a unitcell = atomic volume
            const double atomic_volume = cs::unit_cell.dimensions[0]*cs::unit_cell.dimensions[1]*cs::unit_cell.dimensions[2]/num_atoms_in_unit_cell;

            int exchange_type = 0;

            //what is the exchange type?
            switch(exchange_type){

               case 0: // isotropic

               //loops over all atoms
            //   std::cout << num_local_atoms << std::endl;
               for (int atom = 0; atom <num_local_atoms; atom++){
                  //saves the cell the atom is in and the material
                  const int cell  = cell_array[atom];
                  const int mat   = type_array[atom];
                  //the nearest neighbours are stored in an array - for each atom the start and end index for the array are found,
                  const int start = atoms::neighbour_list_start_index[atom];
                  const int end   = atoms::neighbour_list_end_index[atom] + 1;
              //    std::cout << "cell:\t" << cell << '\t' <<"mat:\t"  << mat << "\t" << start << '\t' << end << std::endl;
                  //loops over all nearest neighbours
                  for(int nn=start;nn<end;nn++){
                     //calcualted the atom id and cell id of the nn atom
                     const int natom = atoms::neighbour_list_array[nn];
                     const int ncell = cell_array[natom];

                     //if interaction is accross a cell boundary

                     if(ncell !=cell){
                        //cacualte the distance between the two atoms
                        const double dx = x_coord_array[atom] - x_coord_array[natom];
                        const double dy = y_coord_array[atom] - y_coord_array[natom];
                        const double dz = z_coord_array[atom] - z_coord_array[natom];
                        const double d2 = dx*dx + dy*dy + dz*dz;

                        //Jij is stored as Jij/mu_s so to get Jij we have to multiply by mu_s
                        //Jij = sum(Jij*distance)
                      //  std::cout << Jij << std::endl;
                        const double Jij=atoms::i_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij*mp::material[mat].mu_s_SI;
                        a2d[cell][ncell] += Jij*d2;
                        num_interactions_cell[cell][ncell] ++;

                     }
                  }
               }

               break;


               case 1: // vector
               terminaltextcolor(RED);
               std::cerr << "Error! Vectoral exchange calculation not yet implemented in micromagnetic mode" << std::endl;
               terminaltextcolor(WHITE);
               zlog << zTs() << "Error! Vectoral exchange calculation not yet implemented in micromagnetic mode" << std::endl;
               err::vexit();
               break;

               case 2: // tensor
               terminaltextcolor(RED);
               std::cerr << "Error! Tensor exchange calculation not yet implemented in micromagnetic mode" << std::endl;
               terminaltextcolor(WHITE);
               zlog << zTs() << "Error! Tensor exchange calculation not yet implemented in micromagnetic mode" << std::endl;
               err::vexit();
               break;
            }

            //unrolls the 2D array of how many interactions accross each boundary into a 1D array.
            int count =0;
            for (int i = 0; i < num_cells; i ++){
               for (int j = 0; j < num_cells; j++){
                  a_test[count] = a2d[i][j];
                  N[count] = num_interactions_cell[i][j];
            //    if (vmpi::my_rank == 0)  std::cerr << i << '\t' << j << '\t' << a2d[i][j] << '\t' << num_interactions_cell[i][j] <<std::endl;
                  count++;
               }
            }

            //Sums the number of interactions on each processor.
            #ifdef MPICF
               MPI_Allreduce(MPI_IN_PLACE, &a_test[0],     num_cells*num_cells,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
               MPI_Allreduce(MPI_IN_PLACE, &N[0],     num_cells*num_cells,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
            #endif

            int i = 0;
            int j = 0;

            //puts the 1D array back into a 2D array
            for (int count = 0; count < num_cells*num_cells; count ++){
               a2d[i][j] = a_test[count];
               num_interactions_cell[i][j] = N[count];
               j++;
               //same!
            //   std::cout << i << '\t' << j << '\t' << a_test[count] << '\t' << num_interactions_cell[i][j] <<std::endl;
               if (j == num_cells) {
                  i ++;
                  j = 0;
               }
            }


            //sums over all interactions to check interaction between cell i j = interaction cell j i
            //non symetric interactions not realistic
            for (int i = 0; i < num_cells; i ++){

               int celli = i;
               for (int j = 0; j < num_cells; j++){
                  int cellj = j;

                  if (int(a2d[celli][cellj]) != int(a2d[cellj][celli])) std::cout << "Error! Non symmetric exchange" <<"\t"  <<  celli << '\t' << cellj << "\t"  <<  a2d[celli][cellj]<<"\t"  <<  a2d[cellj][celli] <<std::endl;
               }
            }

            // loops over all cells to turn the 2D array into a 1D array
            // multiplys A by cell size/2Ms*V_Atomic to ad din the terms of H_Ex
            //removes all the zero interactions by using neighbourlists.
            //The neighbourlists store every interaction as a list. The section of list relevent to each cell is pointed out using the start index and end index.
            if (num_cells > 1){
               for (int celli =0; celli < num_cells; celli++){
                  double cell_size = pow(volume_array[celli],1./3.);                                        //calcualte the size of each cell
                  macro_neighbour_list_start_index[celli] = array_index;                                    //saves the start index for each cell to an array for easy access later
                  double N = volume_array[celli]/atomic_volume;
                  for (int j =0; j <num_cells; j++){
                     int cellj = j;//local_cell_array[j];
                     if (a2d[celli][cellj] != 0 && celli != cellj){
                        macro_neighbour_list_array.push_back(cellj);                                        //if the interaction is non zero add the cell to the neighbourlist
                        a.push_back(-(a2d[celli][cellj]/(4*atomic_volume)));//*cell_size)/(2.0*ms[celli]*atomic_volume*num_interactions_cell[celli][cellj]));
                        //calcualtes the exchange interaction for the cells.                           //the end index is updated for each cell so is given the value for the last cell.
                        
                        a[array_index] = (a[array_index]*2*cell_size)/(ms[celli]*N); //*N instead od numinteraction                        a[array_index] = (a[array_index]*2*cell_size)/(ms[celli]*N); //*N instead od numinteraction
                       // std::cout << celli << '\t' << cellj << '\t' << a[array_index] << '\t' << cell_size << '\t' << ms[celli] << '\t' << N << "\t" << (a[array_index]*2*cell_size)/(ms[celli]*N) << std::endl;
                        macro_neighbour_list_end_index[celli] = array_index;
                        //the end index is updated for each cell so is given the value for the last cell.
                        array_index ++;
                     }
                  }
               }
            }
//            std::cin.get();
            return a;        //returns a 1D vector of the cellular exchange interactions,
         }
      }
   }
