//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2016, Jack B Collings 2021. All rights reserved.
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

// unitcell module headers
#include "internal.hpp"

namespace unitcell{

   //----------------------------------------------------------------------------
   // Function to initialize unitcell module
   //----------------------------------------------------------------------------
   void initialise(unit_cell_t & unit_cell){

   	// check for read-in of unit cell
   	std::string blank="";
   	if(uc::internal::unit_cell_filename.c_str()!=blank){
   		uc::internal::read_unit_cell(unit_cell, uc::internal::unit_cell_filename);
   		return;
   	}

      //------------------------------------------------------------------------
      // Determine which unit cell to build
      //------------------------------------------------------------------------
      if(uc::internal::crystal_structure      == "sc"             ) uc::internal::build_simple_cubic(           unit_cell );
      else if(uc::internal::crystal_structure == "bcc"            ) uc::internal::build_body_centred_cubic(     unit_cell );
      else if(uc::internal::crystal_structure == "bcc-110"        ) uc::internal::build_body_centred_cubic_110( unit_cell );
      else if(uc::internal::crystal_structure == "fcc"            ) uc::internal::build_face_centred_cubic(     unit_cell );
      else if(uc::internal::crystal_structure == "fcc-111"        ) uc::internal::build_face_centred_cubic_111( unit_cell );
      else if(uc::internal::crystal_structure == "hcp"            ) uc::internal::build_hexagonal_close_packed( unit_cell );
      else if(uc::internal::crystal_structure == "heusler"        ) uc::internal::build_heusler(                unit_cell );
      else if(uc::internal::crystal_structure == "honeycomb"      ) uc::internal::build_honeycomb(              unit_cell );
      else if(uc::internal::crystal_structure == "alpha-honeycomb") uc::internal::build_honeycomb_alpha(        unit_cell );
      else if(uc::internal::crystal_structure == "beta-honeycomb")  uc::internal::build_honeycomb_beta(         unit_cell );
      else if(uc::internal::crystal_structure == "kagome"         ) uc::internal::build_kagome(                 unit_cell );
      else if(uc::internal::crystal_structure == "mn2au"          ) uc::internal::build_mn2au(                  unit_cell );
      else if(uc::internal::crystal_structure == "NdFeB"          ) uc::internal::build_NdFeB(                  unit_cell );
      else if(uc::internal::crystal_structure == "rocksalt"       ) uc::internal::build_rock_salt(              unit_cell );
      else if(uc::internal::crystal_structure == "spinel"         ) uc::internal::build_spinel(                 unit_cell );
      else if(uc::internal::crystal_structure == "spinel-layered" ) uc::internal::build_spinel_layered(         unit_cell );
      else if(uc::internal::crystal_structure == "SmFeN"          ) uc::internal::build_SmFeN(                  unit_cell );
      // Otherwise print an error to user
      else{
         terminaltextcolor(RED);
         std::cerr << "Error: Unknown crystal_type "<< uc::internal::crystal_structure << " found during unit cell initialisation. Exiting." << std::endl;
         terminaltextcolor(WHITE);
         zlog << zTs() << "Error: Unknown crystal_type "<< uc::internal::crystal_structure << " found during unit cell initialisation. Exiting." << std::endl;
         err::vexit();
      }

      // optionally write generated unit cell file to disk
      //internal::write_unit_cell_file(unit_cell);

      return;

   }

   //------------------------------------------------------------------------
   // simple function to force simple cubic crystal structure
   //------------------------------------------------------------------------
   void set_crystal_structure_to_simple_cubic(){
      uc::internal::crystal_structure = "sc";
      return;
   }

} // end of unitcell namespace
