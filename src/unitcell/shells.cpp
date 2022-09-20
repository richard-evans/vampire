//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2020. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <algorithm>
#include <cmath>
#include <iostream>

// Vampire headers
#include "errors.hpp"
#include "unitcell.hpp"
#include "vio.hpp"

// unitcell module headers
#include "internal.hpp"

namespace unitcell{

// simple struct to store list of range/ID pairs
struct shell_t{
   unsigned int id; // interaction ID
   unsigned int shell; // shell number
   double range; // range in cell size
};

// comparison function
bool compare(shell_t first, shell_t second){
	if( first.range < second.range) return true;
	else return false;
}

//------------------------------------------------------------------------
// Class function to determine interaction shells
//------------------------------------------------------------------------
void unitcell::exchange_template_t::find_shells(){

   // if shells not needed then do nothing
   if(uc::internal::exchange_function != internal::shell) return;

   // store list of IDs and ranges for sorting
   std::vector<shell_t> interaction_list;

   // load in computed interactions
   for(size_t i=0; i<interaction.size(); i++){
      shell_t tmp;
      tmp.id = i;
      tmp.range = interaction[i].rij;
      interaction_list.push_back(tmp);
   }

   // sort list in range order
   std::sort(interaction_list.begin(), interaction_list.end(), compare);



   //-------------------------------------------
   // determine sets of shells within tolerance
   //-------------------------------------------
   unsigned int shell = 0; // initial shell
   const double tolerance = 0.001; // fractions of unit cell
   double current_range = 0.0;
   if(interaction_list.size() > 0) current_range = interaction_list[0].range; // updating value of shell range
   else{
      std::cerr << "Programmer error!: Interaction list contains no atoms causing seg fault in shell calculation!" << std::endl;
      err::vexit();
   }

   std::vector<int> shell_count(1,0); // list of number of neighbours in each shell
   std::vector<double> shell_range(1,interaction_list[0].range); // list of number of neighbours in each shell

   // check atom at roughly the same range, and if so lump into the same shell
   for(size_t i=0; i<interaction.size(); i++){
      if(interaction_list[i].range < current_range + tolerance){
         interaction_list[i].shell = shell;
         shell_count[shell]++; // increment shell counter
      }
      else{
         // neighbour farther than tolerance and must be in the next shell
         shell++;
         interaction_list[i].shell = shell;
         current_range = interaction_list[i].range;
         // increment records of shell range and counters
         shell_count.push_back(1);
         shell_range.push_back(interaction_list[i].range);
      }
   }

   // Save shell numbers in interaction list
   for(size_t i=0; i<interaction_list.size(); i++){
      int id = interaction_list[i].id;
      interaction[id].shell = interaction_list[i].shell;
   }
   // print
   //for(int i=0; i<interaction.size(); i++){
   //   std::cout << i << "\t" << interaction_list[i].id << "\t" << interaction_list[i].range << "\t" << interaction_list[i].shell << std::endl;
   //}

   // print computed shell information to log file
   zlog << zTs() << "Using long-ranged exchange with computed interaction shells" << std::endl;
   zlog << zTs() << "Number of calculated interaction shells: " << shell_count.size() << std::endl;
   zlog << zTs() << "   Shell \tNumber \tRange \t Cumulative" << std::endl;
   int cumulative = 0;
   const int num_atoms = num_unit_cell_atoms;
   for(size_t i=0; i < shell_count.size(); i++){
      cumulative += shell_count[i];
      zlog << zTs() << "     " << i+1 << "   \t" << shell_count[i]/num_atoms << "\t" << shell_range[i] << " \t" << cumulative/num_atoms << std::endl;
   }

   return;

}

} // end if namespace unitcell
