//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "atoms.hpp"
#include "errors.hpp"
#include "exchange.hpp"
#include "vio.hpp"

// exchange module headers
#include "internal.hpp"

namespace exchange{

   //----------------------------------------------------------------------------
   // Function to initialize exchange module
   //----------------------------------------------------------------------------
   void initialize(std::vector<std::vector <cs::neighbour_t> >& cneighbourlist){

      zlog << zTs() << "Initialising data structures for exchange calculation." << std::endl;

      //-------------------------------------------------
   	//	Calculate total number of neighbours
   	//-------------------------------------------------
   	unsigned int counter = 0;

   	for(int atom=0; atom < atoms::num_atoms; atom++){
   		counter+=cneighbourlist[atom].size();
   	}

   	atoms::total_num_neighbours = counter;

   	atoms::neighbour_list_array.resize(atoms::total_num_neighbours,0);
   	atoms::neighbour_interaction_type_array.resize(atoms::total_num_neighbours,0);
      //atoms::neighbour_eij_array.resize(atoms::total_num_neighbours);

   	atoms::neighbour_list_start_index.resize(atoms::num_atoms,0);
   	atoms::neighbour_list_end_index.resize(atoms::num_atoms,0);

   	//	Populate 1D neighbourlist and index arrays
   	counter = 0;
   	for(int atom=0; atom < atoms::num_atoms; atom++){
   		//std::cout << atom << ": ";
   		// Set start index
   		atoms::neighbour_list_start_index[atom]=counter;
   		for(unsigned int nn=0;nn<cneighbourlist[atom].size();nn++){

            // save atom number to 1D interaction list
   			atoms::neighbour_list_array[counter] = cneighbourlist[atom][nn].nn;

   			if(cneighbourlist[atom][nn].nn > atoms::num_atoms){
   				terminaltextcolor(RED);
   				std::cerr << "Fatal Error - neighbour " << cneighbourlist[atom][nn].nn <<" is out of valid range 0-"
   				<< atoms::num_atoms << " on rank " << vmpi::my_rank << std::endl;
   				//std::cerr << "Atom " << atom << " of MPI type " << catom_array[atom].mpi_type << std::endl;
   				terminaltextcolor(WHITE);
   				err::vexit();
   			}

            // save interaction type to 1D array
   			atoms::neighbour_interaction_type_array[counter] = cneighbourlist[atom][nn].i;

   			//std::cout << cneighbourlist[atom][nn] << " ";
   			counter++;
   		}
   		//std::cout << std::endl;
   		// Set end index
   		atoms::neighbour_list_end_index[atom]=counter-1;
   	}

      // Unroll exchange interactions
      exchange::internal::unroll_exchange_interactions();

      // initialise biquadratic_exchange
      exchange::internal::initialize_biquadratic_exchange();

      // Calculate Dzyaloshinskii-Moriya interactions (must be done after exchange unrolling)
      exchange::internal::calculate_dmi(cneighbourlist);

      return;

   }

} // end of exchange namespace
