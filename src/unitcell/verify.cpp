//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "errors.hpp"
#include "unitcell.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// unitcell module headers
#include "internal.hpp"

namespace unitcell{
namespace internal{

//-------------------------------------------------------------------
//
//   Function to verify symmetry of exchange interactions i->j->i
//
//   A non-symmetric interaction list will not work with the MPI
//   parallelization and makes no physial sense.
//
//-------------------------------------------------------------------
void verify_exchange_interactions(unit_cell_t & unit_cell, std::string filename){

   // list of assymetric interactions
   std::vector<int> asym_interaction_list(0);

   // Parallelise in case of large interaction sizes
   int my_num_interactions = unit_cell.interaction.size()/vmpi::num_processors;
   int first = vmpi::my_rank * my_num_interactions;
   int last = first + my_num_interactions;
   if(vmpi::my_rank == vmpi::num_processors-1) last = unit_cell.interaction.size(); // add last points to last processor

   // loop over all interactions to find matching reciprocal interaction
   for(unsigned int i = first; i < last; ++i){

      // Output progress indicator to screen for large interaction counts
      if( (i % (my_num_interactions/10)) == 0 && unit_cell.interaction.size() > 10000) std::cout << "." << std::flush;

      // calculate reciprocal interaction
      unsigned int ia = unit_cell.interaction[i].j;
      unsigned int ja = unit_cell.interaction[i].i;
      int dx = -unit_cell.interaction[i].dx;
      int dy = -unit_cell.interaction[i].dy;
      int dz = -unit_cell.interaction[i].dz;

      // set flag to test for match
      bool match=false;

      // loop over all interactions for reciprocal interactions i -> j -> i
      for(unsigned int j=0; j<unit_cell.interaction.size(); ++j){
         if(unit_cell.interaction[j].i==ia && unit_cell.interaction[j].j==ja && unit_cell.interaction[j].dx==dx && unit_cell.interaction[j].dy==dy && unit_cell.interaction[j].dz==dz){
            match=true;
            break;
         }
      }

      // if no match is found add to list of assymetric interactions
      if(!match){
         asym_interaction_list.push_back(i);
      }
   }

   // Make all processors wait here
   vmpi::barrier();

   // Output error message and list of interactions if found on any processor
   if(asym_interaction_list.size()>0){
      terminaltextcolor(RED);
      std::cerr << "Error! Exchange interaction list in unit cell file " << filename << " contains the following assymetric interactions:" << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << "Error! Exchange interaction list in unit cell file " << filename << " contains the following assymetric interactions:" << std::endl;
      for(unsigned int i=0; i < asym_interaction_list.size(); ++i){
         int id=asym_interaction_list[i];
         terminaltextcolor(RED);
         std::cerr << id << "\t" << unit_cell.interaction[id].i << "\t" << unit_cell.interaction[id].j << "\t" << unit_cell.interaction[id].dx << "\t" << unit_cell.interaction[id].dy << "\t" << unit_cell.interaction[id].dz << std::endl;
         terminaltextcolor(WHITE);
         zlog << "\t\t\t" << id << "\t" << unit_cell.interaction[id].i << "\t" << unit_cell.interaction[id].j << "\t" << unit_cell.interaction[id].dx << "\t" << unit_cell.interaction[id].dy << "\t" << unit_cell.interaction[id].dz << std::endl;
      }
      terminaltextcolor(RED);
      std::cerr << "Assymetric interactions are unphysical: please fix the unit cell file ensuring all interactions are symmetric. Exiting." << std::endl;
      terminaltextcolor(WHITE);
      zlog << "\t\t\t" << "Assymetric interactions are unphysical: please fix the unit cell file ensuring all interactions are symmetric. Exiting." << std::endl;
      err::vexit();
   }

   return;

}

} // end of internal namespace
} // end of unitcell namespace
