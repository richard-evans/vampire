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
#include "atoms.hpp"
#include "errors.hpp"
#include "vio.hpp"

namespace micromagnetic {
namespace internal {

//------------------------------------------------------------------------------
// calculates the curie temperature of each cell
//------------------------------------------------------------------------------
std::vector<double> calculate_tc(int num_local_cells,
                                 std::vector <int> local_cell_array,
                                 int num_atoms,
                                 int num_cells,
                                 std::vector<int> cell_array,                   //1D array storing which cell each atom is in
                                 std::vector<int> neighbour_list_array,         //1D list of the interactions between atoms,
                                 std::vector<int> neighbour_list_start_index,   //the list for each atom in neighbour list array starts at the start index
                                 std::vector<int> neighbour_list_end_index,     //and ends at the end index - stored in these vectors
                                 const std::vector<int> type_array,            //1D array storing which material each atom is
                                 std::vector <mp::materials_t> material){      //class of material parameters for the atoms

   std::vector<double>  J(num_cells,0.0);         //stores the exchange constant for each cell
   std::vector<double>  N(num_cells,0.0);         //1D vector containg the number of atoms in each cell.
   std::vector<double>  Tc(num_cells,0.0);        //1D vector storing the curie temperature of each cell

   const double e = 0.8;                          // mean field constant
   const double kB = 1.38064852e-23;              // boltzman constant

   #ifdef MPICF
       num_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
   #endif

   //-----------------------------------------------------------------------
   //             TC = sum(Jij)* e/(3kB N)
   //------------------------------------------------------------------------

   int exchange_type = 0; //temporary stand-in for exchange type

   switch(exchange_type){

 		case 0: // isotropic
         for (int atom = 0; atom < num_atoms; atom++){
            const int cell  = cell_array[atom];
            const int start = atoms::neighbour_list_start_index[atom];
            const int end   = atoms::neighbour_list_end_index[atom] + 1;
            const int mat   = type_array[atom];
            N[cell]++;
            for(int nn=start;nn<end;nn++){
               const double Jij=atoms::i_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij*mp::material[mat].mu_s_SI;
               J[cell] = J[cell] + Jij;
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

   // Reduce sum of Jij and N on all processors
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &J[0], num_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &N[0], num_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   #endif

   // Set Tc value for all cells
   for (int cell = 0; cell < num_cells; cell++){
      Tc[cell] = J[cell]*e/(3*kB*N[cell]);
   }

   // Reduce Tc for all cells on all processors
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &Tc[0], num_cells, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   #endif

   return Tc;             //returns a 1D array containing the curie temepratures
}

} //closes the internal namspace
}  //closes the micromagnetic namespace
