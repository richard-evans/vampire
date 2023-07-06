//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2018. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "errors.hpp"
#include "material.hpp"
#include "unitcell.hpp"
#include "vio.hpp"

// unitcell module headers
#include "internal.hpp"

namespace unitcell{

   void unitcell::exchange_template_t::read_interactions(
      //unit_cell_t& unit_cell,
      const int num_atoms, // num atoms in unit cell
      std::stringstream& ucf_file,
      std::istringstream& ucf_ss,
      std::string& filename,
      unsigned int& line_counter,
      unsigned int& interaction_range
   ){

		int num_interactions = 0; // assume no interactions
      std::string exchange_type_string; // string defining exchange type

      // get number of exchange types
      ucf_ss >> num_interactions >> exchange_type_string;

      // process exchange string to set exchange type and normalisation
      const int num_exchange_values = unitcell::exchange_template_t::set_exchange_type(exchange_type_string);

      if(num_interactions>=0) interaction.resize(num_interactions);
      else {
         terminaltextcolor(RED);
         std::cerr << "Error! Requested number of biquadratic exchange interactions " << num_interactions << " on line " << line_counter
         << " of unit cell input file " << filename.c_str() << " is less than 0. Exiting" << std::endl; err::vexit();
         terminaltextcolor(WHITE);
      }

      std::cout << "done!\nProcessing data from " << num_interactions << " interactions..." << std::flush;
      zlog << zTs() << "\t" << "Processing unit cell atoms completed" << std::endl;
      zlog << zTs() << "\t" << "Processing data from " << num_interactions << " interactions..." << std::endl;

      // resize interaction counter array to store number of interactions per atom
      ni.resize(num_atoms, 0);

      // loop over all interactions and read into class
      for (int i=0; i<num_interactions; i++){

         // Output progress counter to screen for large interaction counts
         if( (i % (num_interactions/10 + 1)) == 0 && num_interactions > 10000) std::cout << "." << std::flush;

         // declare safe temporaries for interaction input
         int id=i;
         int iatom=-1,jatom=-1; // atom pairs
         int dx=0, dy=0,dz=0; // relative unit cell coordinates
         // get line
         std::string int_line;
         getline(ucf_file,int_line);
         //std::cout << int_line.c_str() << std::endl;
         std::istringstream int_iss(int_line,std::istringstream::in);
         int_iss >> id >> iatom >> jatom >> dx >> dy >> dz;
         //inputfile >> id >> iatom >> jatom >> dx >> dy >> dz;
         line_counter++;
         // check for sane input
         if(iatom>=0 && iatom < num_atoms) interaction[i].i=iatom;
         else if(iatom>=0 && iatom >= num_atoms){
            terminaltextcolor(RED);
            std::cerr << std::endl << "Error! iatom number "<< iatom <<" for interaction id " << id << " on line " << line_counter
                 << " of unit cell input file " << filename.c_str() << " is outside of valid range 0-"
                 << num_atoms-1 << ". Exiting" << std::endl;
            terminaltextcolor(WHITE);
            zlog << zTs() << "Error! iatom number "<< iatom <<" for interaction id " << id << " on line " << line_counter
                 << " of unit cell input file " << filename.c_str() << " is outside of valid range 0-"<< num_atoms-1
                 << ". Exiting" << std::endl;
            err::vexit();
         }
         else{
           terminaltextcolor(RED);
           std::cerr << std::endl << "Error! No valid interaction for interaction id " << id << " on line " << line_counter
                << " of unit cell input file " << filename.c_str() << ". Possibly too many interactions defined. Exiting" << std::endl;
           terminaltextcolor(WHITE);
           zlog << zTs() << "Error! No valid interaction for interaction id " << id << " on line " << line_counter
                << " of unit cell input file " << filename.c_str() << ". Possibly too many interactions defined. Exiting" << std::endl;
           err::vexit();
         }
         if(iatom>=0 && jatom < num_atoms) interaction[i].j=jatom;
         else{
            terminaltextcolor(RED);
            std::cerr << std::endl << "Error! jatom number "<< jatom <<" for interaction id " << id << " on line " << line_counter
                 << " of unit cell input file " << filename.c_str() << " is outside of valid range 0-"
                 << num_atoms-1 << ". Exiting" << std::endl;
            terminaltextcolor(WHITE);
            zlog << zTs() << "Error! jatom number "<< jatom <<" for interaction id " << id << " on line " << line_counter
                 << " of unit cell input file " << filename.c_str() << " is outside of valid range 0-"
                 << num_atoms-1 << ". Exiting" << std::endl;
            err::vexit();
         }
         interaction[i].dx=dx;
         interaction[i].dy=dy;
         interaction[i].dz=dz;

         // check for long range interactions
         if(static_cast<unsigned int>(abs(dx))>interaction_range) interaction_range=abs(dx);
         if(static_cast<unsigned int>(abs(dy))>interaction_range) interaction_range=abs(dy);
         if(static_cast<unsigned int>(abs(dz))>interaction_range) interaction_range=abs(dz);

         //int iatom_mat = unit_cell.atom[iatom].mat;
         //int jatom_mat = unit_cell.atom[jatom].mat;
         switch(num_exchange_values){
            case 1:
               int_iss >> interaction[i].Jij[0][0];
               // save interactions into diagonal components of the exchange tensor
               interaction[i].Jij[1][1] = interaction[i].Jij[0][0];
               interaction[i].Jij[2][2] = interaction[i].Jij[0][0];
               break;
            case 3:
               int_iss >> interaction[i].Jij[0][0] >> interaction[i].Jij[1][1] >> interaction[i].Jij[2][2];
               break;
            case 9:
               int_iss >> interaction[i].Jij[0][0] >> interaction[i].Jij[0][1] >> interaction[i].Jij[0][2];
               int_iss >> interaction[i].Jij[1][0] >> interaction[i].Jij[1][1] >> interaction[i].Jij[1][2];
               int_iss >> interaction[i].Jij[2][0] >> interaction[i].Jij[2][1] >> interaction[i].Jij[2][2];
               break;
            default:
               terminaltextcolor(RED);
               std::cerr << "Programmer Error! Requested number of exchange values " << num_exchange_values << " on line " << line_counter
         << " of unit cell input file " << filename.c_str() << " is outside of valid range 1,3 or 9. Exiting" << std::endl;
               terminaltextcolor(WHITE);
               err::vexit();
         }

         // increment number of interactions for atom i
         ni[iatom]++;

      }

      return;

   }

} // end of unitcell namespace
