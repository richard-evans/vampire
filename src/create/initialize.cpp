//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <iostream>

// Vampire headers
#include "create.hpp"
#include "internal.hpp"
#include "vio.hpp"

namespace create{

   //-------------------------------------------------------------------------------
   // Function to initialize create module
   //-------------------------------------------------------------------------------
   void initialize(){

      // If create material parameters uninitialised then initialise with default parameters
      if(create::internal::mp.size() == 0){
         create::internal::mp.resize(mp::num_materials);
      }

      // Loop over materials to check for invalid input and warn appropriately
		for(int mat=0;mat<mp::num_materials;mat++){
			const double lmin=create::internal::mp[mat].min;
			const double lmax=create::internal::mp[mat].max;
			for(int nmat=0;nmat<mp::num_materials;nmat++){
				if(nmat!=mat){
					double min=create::internal::mp[nmat].min;
					double max=create::internal::mp[nmat].max;
					if(((lmin>min) && (lmin<max)) || ((lmax>min) && (lmax<max))){
						terminaltextcolor(RED);
						std::cout << "Warning: Overlapping material heights found. Check log for details." << std::endl;
						terminaltextcolor(WHITE);
						zlog << zTs() << "Warning: material " << mat+1 << " overlaps material " << nmat+1 << "." << std::endl;
						zlog << zTs() << "If you have defined geometry then this may be OK, or possibly you meant to specify alloy keyword instead." << std::endl;
						zlog << zTs() << "----------------------------------------------------" << std::endl;
						zlog << zTs() << "  Material "<< mat+1 << ":minimum-height = " << lmin << std::endl;
						zlog << zTs() << "  Material "<< mat+1 << ":maximum-height = " << lmax << std::endl;
						zlog << zTs() << "  Material "<< nmat+1 << ":minimum-height = " << min << std::endl;
						zlog << zTs() << "  Material "<< nmat+1 << ":maximum-height = " << max << std::endl;
					}
				}
			}
		}
      for (int mat = 0; mat < mp::num_materials; mat++){
         mp::material[mat].pinning_field_unit_vector[0] *= mp::material[mat].pinning_field_strength;
         mp::material[mat].pinning_field_unit_vector[1] *= mp::material[mat].pinning_field_strength;
         mp::material[mat].pinning_field_unit_vector[2] *= mp::material[mat].pinning_field_strength;
      }

      //------------------------------------------------------------------------
      // output informative message on periodic boundaries to log file
      //------------------------------------------------------------------------
      std::string xstr = cs::pbc[0] ? (cs::pbc[1] || cs::pbc[2] ? "x," : "x") : "";
      std::string ystr = cs::pbc[1] ? (cs::pbc[2]               ? "y," : "y") : "";
      std::string zstr = cs::pbc[2] ? "z" : "";
      // set plural only if at least two directions are set
      std::string plural = (cs::pbc[0] && cs::pbc[1]) || (cs::pbc[0] && cs::pbc[2]) || (cs::pbc[1] && cs::pbc[2]) ? "s" : "";
      // only if sirections set output message to user
      if(cs::pbc[0] || cs::pbc[1] || cs::pbc[2] ) zlog << zTs() << "Setting periodic boundary conditions in " << xstr << ystr << zstr << " direction" << plural << std::endl;



      return;
   }

} // end of namespace create
