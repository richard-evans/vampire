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
   void initialize(std::vector<std::vector <neighbours::neighbour_t> >& bilinear,
                   std::vector<std::vector <neighbours::neighbour_t> >& biquadratic){

      zlog << zTs() << "Initialising data structures for exchange calculation." << std::endl;

      //-------------------------------------------------
   	//	Calculate total number of neighbours
   	//-------------------------------------------------
   	uint64_t counter = 0; // number of exchange interactions

   	for(uint64_t atom = 0; atom < atoms::num_atoms; atom++){
   		counter += bilinear[atom].size();
   	}

   	atoms::total_num_neighbours = counter;

   	atoms::neighbour_list_array.resize(atoms::total_num_neighbours,0);
   	atoms::neighbour_interaction_type_array.resize(atoms::total_num_neighbours,0);
      //atoms::neighbour_eij_array.resize(atoms::total_num_neighbours);

   	atoms::neighbour_list_start_index.resize(atoms::num_atoms,0);
   	atoms::neighbour_list_end_index.resize(atoms::num_atoms,0);

   	//	Populate 1D neighbourlist and index arrays
   	counter = 0;
   	for(uint64_t atom=0; atom < atoms::num_atoms; atom++){
   		//std::cout << atom << ": ";
   		// Set start index
   		atoms::neighbour_list_start_index[atom]=counter;
   		for(unsigned int nn=0;nn<bilinear[atom].size();nn++){

            // save atom number to 1D interaction list
   			atoms::neighbour_list_array[counter] = bilinear[atom][nn].nn;

   			if(bilinear[atom][nn].nn > atoms::num_atoms){
   				terminaltextcolor(RED);
   				std::cerr << "Fatal Error - neighbour " << bilinear[atom][nn].nn <<" is out of valid range 0-"
   				<< atoms::num_atoms << " on rank " << vmpi::my_rank << std::endl;
   				//std::cerr << "Atom " << atom << " of MPI type " << catom_array[atom].mpi_type << std::endl;
   				terminaltextcolor(WHITE);
   				err::vexit();
   			}

            // save interaction type to 1D array
   			atoms::neighbour_interaction_type_array[counter] = bilinear[atom][nn].i;

   			//std::cout << bilinear[atom][nn] << " ";
   			counter++;
   		}
   		//std::cout << std::endl;
   		// Set end index
   		atoms::neighbour_list_end_index[atom]=counter-1;
   	}

      //------------------------------------------------------------------------
      // Biquadratic exchange list - should probably be classified at some point
      //------------------------------------------------------------------------
      if(exchange::biquadratic){
         counter = 0;
         for(uint64_t atom = 0; atom < atoms::num_atoms; atom++){
      		counter += biquadratic[atom].size();
      	}

      	exchange::internal::biquadratic_neighbour_list_array.resize(counter,0);
      	exchange::internal::biquadratic_neighbour_interaction_type_array.resize(counter,0);

      	exchange::internal::biquadratic_neighbour_list_start_index.resize(atoms::num_atoms,0);
      	exchange::internal::biquadratic_neighbour_list_end_index.resize(atoms::num_atoms,0);

      	//	Populate 1D neighbourlist and index arrays
      	counter = 0;
      	for(int atom=0; atom < atoms::num_atoms; atom++){
      		//std::cout << atom << ": ";
      		// Set start index
      		exchange::internal::biquadratic_neighbour_list_start_index[atom]=counter;
      		for(unsigned int nn=0; nn < biquadratic[atom].size(); nn++){

               // save atom number to 1D interaction list
      			exchange::internal::biquadratic_neighbour_list_array[counter] = biquadratic[atom][nn].nn;

      			if(biquadratic[atom][nn].nn > atoms::num_atoms){
      				terminaltextcolor(RED);
      				std::cerr << "Fatal Error - biquadratic neighbour " << biquadratic[atom][nn].nn <<" is out of valid range 0-"
      				<< atoms::num_atoms << " on rank " << vmpi::my_rank << std::endl;
      				//std::cerr << "Atom " << atom << " of MPI type " << catom_array[atom].mpi_type << std::endl;
      				terminaltextcolor(WHITE);
      				err::vexit();
      			}

               // save interaction type to 1D array
      			exchange::internal::biquadratic_neighbour_interaction_type_array[counter] = biquadratic[atom][nn].i;

      			//std::cout << biquadratic[atom][nn] << " ";
      			counter++;

      		}
      		//std::cout << std::endl;
      		// Set end index
      		exchange::internal::biquadratic_neighbour_list_end_index[atom]=counter-1;
      	}
      }

      // Unroll exchange interactions
      exchange::internal::unroll_exchange_interactions();

      // initialise biquadratic_exchange
      exchange::internal::initialize_biquadratic_exchange();

      // Calculate Dzyaloshinskii-Moriya interactions (must be done after exchange unrolling)
      exchange::internal::calculate_dmi(bilinear);

      return;

   }

} // end of exchange namespace
