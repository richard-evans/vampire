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

   // Error function when non-symmetric exchange interactions are detected
   void exchange_error(int i, int j, int s, std::vector<double> Jij, std::vector<double> Jji, std::string descriptor){

      // Error found - report to user (in terms of material ID and shell and terminate program
      terminaltextcolor(RED);
         std::cerr << "Error! Non-symmetric " << descriptor << " exchange interactions for materials " << i+1 << " and " << j+1 << " in shell " << s+1 << "." << std::endl;
         std::cerr << "See log file for details. Exiting ungracefully." << std::endl;
      terminaltextcolor(WHITE);

      // report detailed info in the log file
      zlog << zTs() << "Error! Non-symmetric " << descriptor << "exchange interactions for materials " << i+1 << " and " << j+1 << " in shell " << s+1 << std::endl;
      zlog << zTs() << "\tmaterial[" << i+1 << "]:" << descriptor << "-exchange-matrix-" << s+1 << "-nn[" << j+1 << "] = ";
         for(size_t vi = 0; vi < Jij.size(); vi++) zlog << Jij[vi] << "\t";
         zlog << std::endl;
      zlog << zTs() << "\tmaterial[" << j+1 << "]:" << descriptor << "-exchange-matrix-" << s+1 << "-nn[" << i+1 << "] = ";
         for(size_t vi = 0; vi < Jji.size(); vi++) zlog << Jji[vi] << "\t";
         zlog << std::endl;
      zlog << zTs() << "\tThe definition of " << descriptor << " exchange requires that there are the same number of values and that these values are the same. Exiting." << std::endl;

      // exit program ungracefully
      err::vexit();

   }

   //----------------------------------------------------------------------------
   // Function to initialize exchange module
   //----------------------------------------------------------------------------
   void initialize(std::vector<std::vector <neighbours::neighbour_t> >& bilinear,
                   std::vector<std::vector <neighbours::neighbour_t> >& biquadratic){

      zlog << zTs() << "Initialising data structures for exchange calculation." << std::endl;

      //-----------------------------------------------------------------------------
      // Optionally output exchange template and halo
      //-----------------------------------------------------------------------------

      //std::cout << "Num atoms: " << atoms::num_atoms << std::endl;
      //std::cout << "Num core atoms: " << vmpi::num_core_atoms << std::endl;
      //std::cout << "Num boundary atoms: " << vmpi::num_bdry_atoms << std::endl;
      //std::cout << "Num halo atoms: " << vmpi::num_halo_atoms << std::endl;

      //std::ofstream ofile("mpi_crystal_after.xyz");

      //ofile << atoms::num_atoms << "\n\n";
      //for(uint64_t atom=0; atom < atoms::num_atoms; atom++){
      //   std::string ele = "Ag";
      //   if(atom < vmpi::num_core_atoms) ele = "Ag";
      //   else if(atom < vmpi::num_core_atoms+vmpi::num_bdry_atoms) ele = "Fe";
      //   else if(atom < vmpi::num_core_atoms + vmpi::num_bdry_atoms + vmpi::num_halo_atoms) ele = "Li";
      //   ofile << ele << "\t" << atoms::x_coord_array[atom] << "\t" << atoms::y_coord_array[atom] << "\t" << atoms::z_coord_array[atom] << std::endl;
      //}

      //ofile.close();

      //-----------------------------------------------------------------------------
      // Check that bilinear exchange matrix is symmetric for i->j and j->i interactions
      //-----------------------------------------------------------------------------
      const int max_materials = exchange::internal::bilinear_exchange_constants.get_max_materials();
      const int max_shells = exchange::internal::bilinear_exchange_constants.get_max_shells();
      // loop over half of exchange matrix of set values to check symmetry
      for(int i = 0; i < max_materials; i++){
         for(int j = 0; j < i; j++){
            for(int s = 0; s < max_shells; s++){
               // get i->j exchange
               std::vector<double> Jij = exchange::internal::bilinear_exchange_constants.get_exchange_values(i,j,s);
               // get j->i exchange
               std::vector<double> Jji = exchange::internal::bilinear_exchange_constants.get_exchange_values(j,i,s);
               // check for same number of exchange values
               if( Jij.size() != Jji.size() ) exchange_error(i, j, s, Jij, Jji, "");

               // check for rational values for all elements
               for(size_t v = 0 ; v < Jij.size() ; v++){

                  // Check for non-zero value (avoids divide by zero)
                  if(fabs( Jji[v] ) > 0.0){

                     // Calculate ratio of i->j / j-> exchange constants
                     double ratio = Jij[v] / Jji[v];

                     // Check that ratio ~ 1.0 for symmetric exchange interactions
                     if( (ratio < 0.99999) || (ratio > 1.00001) ) exchange_error(i, j , s, Jij, Jji, "");

                  }
                  else{
                     // Jji is zero, now check for Jij == 0.0+tolerance
                     if(fabs( Jij[v] ) > 1.0e-200) exchange_error(i, j, s, Jij, Jji, "");
                  }
               }
            }
         }
      }

      //-----------------------------------------------------------------------------
      // Check that bilinear exchange matrix is symmetric for i->j and j->i interactions
      //-----------------------------------------------------------------------------
      const int max_bq_materials = exchange::internal::biquadratic_exchange_constants.get_max_materials();
      const int max_bq_shells    = exchange::internal::biquadratic_exchange_constants.get_max_shells();
      // loop over half of exchange matrix of set values to check symmetry
      for(int i = 0; i < max_bq_materials; i++){
         for(int j = 0; j < i; j++){
            for(int s = 0; s < max_bq_shells; s++){
               // get i->j exchange
               std::vector<double> Jij = exchange::internal::biquadratic_exchange_constants.get_exchange_values(i,j,s);
               // get j->i exchange
               std::vector<double> Jji = exchange::internal::biquadratic_exchange_constants.get_exchange_values(j,i,s);
               // check for same number of exchange values
               if( Jij.size() != Jji.size() ) exchange_error(i, j, s, Jij, Jji, "biquadratic");

               // check for rational values for all elements
               for(size_t v = 0 ; v < Jij.size() ; v++){

                  // Check for non-zero value (avoids divide by zero)
                  if(fabs( Jji[v] ) > 0.0){

                     // Calculate ratio of i->j / j-> exchange constants
                     double ratio = Jij[v] / Jji[v];

                     // Check that ratio ~ 1.0 for symmetric exchange interactions
                     if( (ratio < 0.99999) || (ratio > 1.00001) ) exchange_error(i, j , s, Jij, Jji, "biquadratic");

                  }
                  else{
                     // Jji is zero, now check for Jij == 0.0+tolerance
                     if(fabs( Jij[v] ) > 1.0e-200) exchange_error(i, j, s, Jij, Jji, "biquadratic");
                  }
               }
            }
         }
      }

      //-----------------------------------------------------------------------------
      // save type of interaction template and if material file constants are used
      exchange::internal::exchange_type = cs::unit_cell.bilinear.exchange_type;
      exchange::internal::use_material_exchange_constants = cs::unit_cell.bilinear.use_material_exchange_constants;

      // check that minimum required exchange type is enabled
      if(exchange::internal::exchange_type == exchange::isotropic &&  exchange::internal::minimum_needed_exchange_type == exchange::vectorial){
         exchange::internal::exchange_type = exchange::vectorial;
      }

      //-------------------------------------------------
   	//	Calculate total number of neighbours
   	//-------------------------------------------------
   	uint64_t counter = 0; // number of exchange interactions

   	for(uint64_t atom = 0; atom < static_cast<uint64_t>(atoms::num_atoms); atom++){
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
   	for(uint64_t atom=0; atom < static_cast<uint64_t>(atoms::num_atoms); atom++){
   		//std::cout << atom << ": ";
   		// Set start index
   		atoms::neighbour_list_start_index[atom]=counter;
   		for(size_t nn=0;nn<bilinear[atom].size();nn++){

            // save atom number to 1D interaction list
   			atoms::neighbour_list_array[counter] = bilinear[atom][nn].nn;

   			if(bilinear[atom][nn].nn >= atoms::num_atoms){
   				terminaltextcolor(RED);
   				std::cerr << "Fatal Error - neighbour atom " << bilinear[atom][nn].nn <<" is out of valid range 0-"
   				<< atoms::num_atoms-1 << " on rank " << vmpi::my_rank << std::endl;
               std::cerr << "\tAtom number      : " << atom << std::endl;
               std::cerr << "\tNeighbour number : " << nn << std::endl;
               std::cerr << "\tNeighbour atom   : " << atoms::neighbour_list_array[counter] << std::endl;
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

         // determine total number of biquadratic exchange interactions
         counter = 0;
         for(uint64_t atom = 0; atom < static_cast<uint64_t>(atoms::num_atoms); atom++){
      		counter += biquadratic[atom].size();
      	}

         // save type of interaction template and if material file constants are used
         exchange::internal::biquadratic_exchange_type = cs::unit_cell.biquadratic.exchange_type;
         exchange::internal::use_material_biquadratic_exchange_constants = cs::unit_cell.biquadratic.use_material_exchange_constants;

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

      			if(biquadratic[atom][nn].nn >= atoms::num_atoms){
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
      exchange::internal::unroll_exchange_interactions(bilinear);

      // initialise biquadratic_exchange
      exchange::internal::initialize_biquadratic_exchange();

      // Calculate Dzyaloshinskii-Moriya interactions (must be done after exchange unrolling)
      exchange::internal::calculate_dmi(bilinear);

      // Calculate Kitaev interactions (must be done after exchange unrolling)
      exchange::internal::calculate_kitaev(bilinear);

      return;

   }

} // end of exchange namespace
