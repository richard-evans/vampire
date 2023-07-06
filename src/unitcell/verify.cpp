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

//-------------------------------------------------------------------
//
//   Function to verify symmetry of exchange interactions i->j->i
//
//   A non-symmetric interaction list will not work with the MPI
//   parallelization and makes no physical sense.
//
//-------------------------------------------------------------------
void unitcell::exchange_template_t::verify(std::string filename){

   // ignore any cases with more than 100,000 interactions
   if(interaction.size() > 100000) return;

   // list of assymetric interactions
   std::vector<int> asym_interaction_list(0);

   // Parallelise in case of large interaction sizes
   size_t my_num_interactions = interaction.size()/vmpi::num_processors;
   size_t first = vmpi::my_rank * my_num_interactions;
   size_t last = first + my_num_interactions;
   if(vmpi::my_rank == vmpi::num_processors-1) last = interaction.size(); // add last points to last processor

   // loop over all interactions to find matching reciprocal interaction
   for(size_t i = first; i < last; ++i){

      // Output progress indicator to screen for large interaction counts
      if( (i % (my_num_interactions/10 + 1)) == 0 && interaction.size() > 10000) std::cout << "." << std::flush;

      // calculate reciprocal interaction
      unsigned int ia = interaction[i].j;
      unsigned int ja = interaction[i].i;
      int dx = -interaction[i].dx;
      int dy = -interaction[i].dy;
      int dz = -interaction[i].dz;

      // set flag to test for match
      bool match=false;

      // loop over all interactions for reciprocal interactions i -> j -> i
      for(size_t j=0; j<interaction.size(); ++j){
         if(interaction[j].i==ia && interaction[j].j==ja && interaction[j].dx==dx && interaction[j].dy==dy && interaction[j].dz==dz){
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
         std::cerr << id << "\t" << interaction[id].i << "\t" << interaction[id].j << "\t" << interaction[id].dx << "\t" << interaction[id].dy << "\t" << interaction[id].dz << std::endl;
         terminaltextcolor(WHITE);
         zlog << "\t\t\t" << id << "\t" << interaction[id].i << "\t" << interaction[id].j << "\t" << interaction[id].dx << "\t" << interaction[id].dy << "\t" << interaction[id].dz << std::endl;
      }
      terminaltextcolor(RED);
      std::cerr << "Assymetric interactions are unphysical: please fix the unit cell file ensuring all interactions are symmetric. Exiting." << std::endl;
      terminaltextcolor(WHITE);
      zlog << "\t\t\t" << "Assymetric interactions are unphysical: please fix the unit cell file ensuring all interactions are symmetric. Exiting." << std::endl;
      err::vexit();
   }

   return;

}

} // end of unitcell namespace
