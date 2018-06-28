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
#include "errors.hpp"
#include "vio.hpp"

namespace micromagnetic{

   namespace internal{


      //function to calulate the intercell exchange.
      //Each cell has an exchange interaction with every other cell in the system.
      //This is calculated as a sum over all the atoms on the surface of the cells and there interaction with atoms in the other cell.

      //---------------------------------------------------------------------------------------------------
      //
      //                       A = 1/4*sum((Jij) *(x_i - x_j)^2)
      //This calculates the interaction between neighbouring cells therefore only interactions where i and j
      //are in different macrocells are summed.
      //
      // we take some components of the exchange field in here and calulate them as part of a therefore we return
      //
      //      A = (1/2)*sum((Jij)(x_i - x_j)^2) * (V_cell/V_atomic) * (1/l_cell) * (1/Ms)
      //
      //---------------------------------------------------------------------------------------

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

            std::vector<double> a; //1D vector to store the 2D vector odf the exchange intereactions
            a.resize(num_cells*num_cells,0.0);

            //calculates the atomic volume  = volume of one cell/number of atoms in a unitcell = atomic volume
            const double atomic_volume = cs::unit_cell.dimensions[0]*cs::unit_cell.dimensions[1]*cs::unit_cell.dimensions[2]/num_atoms_in_unit_cell;

            int exchange_type = 0;

            //what is the exchange type?
            switch(exchange_type){

               case 0: // isotropic

               //loops over all atoms
               for (int atom = 0; atom <num_atoms; atom++){
                  //saves the cell the atom is in and the material
                  const int cell  = cell_array[atom];
                  const int mat   = type_array[atom];
               //   std::cout << cell << '\t' << mat << std::endl;
                  //the nearest neighbours are stored in an array - for each atom the start and end index for the array are found,
                  const int start = atoms::neighbour_list_start_index[atom];
                  const int end   = atoms::neighbour_list_end_index[atom] + 1;
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
                     //   std::cout << d2 <<std::endl;
                     //   std::cin.get();
                        //Jij is stored as Jij/mu_s so to get Jij we have to multiply by mu_s
                        //Jij = sum(Jij*distance)
                        const double Jij=atoms::i_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij*mp::material[mat].mu_s_SI;
                        a2d[cell][ncell] += Jij*d2;
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

            //sums over all interactions to check interaction between cell i j = interaction cell j i
            //non symetric interactions not realistic
            for (int i = 0; i < num_local_cells; i ++){
               int celli = i;//local_cell_array[i];
               for (int j = 0; j < num_local_cells; j++){
                  int cellj = local_cell_array[i];
                  if (int(a2d[celli][cellj]) != int(a2d[cellj][celli])) std::cout << "Error! Non symetric exchange" <<"\t"  <<  celli << '\t' << cellj << "\t"  <<  a2d[celli][cellj]<<"\t"  <<  a2d[cellj][celli] <<std::endl;
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

               for (int j =0; j <num_local_cells; j++){
                  int cellj = local_cell_array[j];
                  if (a2d[celli][cellj] != 0){
                     macro_neighbour_list_array.push_back(cellj);                                        //if the interaction is non zero add the cell to the neighbourlist
                     a[array_index] = -(a2d[celli][cellj]*cell_size)/(2*ms[celli]*atomic_volume);                 //calcualtes the exchange interaction for the cells.
                     macro_neighbour_list_end_index[celli] = array_index;
                  //   std::cout << celli << '\t' << cellj << '\t' <<   a[array_index] <<std::endl;                           //the end index is updated for each cell so is given the value for the last cell.
                     array_index ++;
                  }
               }


            }
         }
            //   #ifdef MPICF
            //     MPI_Allreduce(MPI_IN_PLACE, &a[0],     array_index,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
            //   #endif

            return a;        //returns a 1D vector of the cellular exchange interactions,
         }
      }
   }
