//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard Evans 2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// vampire headers
#include "create.hpp"
#include "vio.hpp"

// create module headers
#include "internal.hpp"

namespace create{

namespace internal{

//------------------------------------------------------------------------------
// Simple function to calculate atomic composition
//------------------------------------------------------------------------------
void calculate_atomic_composition(std::vector<cs::catom_t> & catom_array){

   // print informative message to log file
	zlog<< zTs() << "Determining atomic composition:" << std::endl;

   // temporary counter for number of (local) atoms
   unsigned int num_atoms = 0;

   // array to store number of atoms in each material
   std::vector<uint64_t> material_numbers(mp::num_materials,0);

   // loop over all atoms to calculate number of atoms in each material
   for(unsigned int atom = 0; atom < catom_array.size(); atom++){

      //std::cout << atom << "\t" << catom_array[atom].mpi_type << std::endl;

      // only count local (core and boundary) atoms
      if(catom_array[atom].mpi_type !=2 ){
         material_numbers.at(catom_array[atom].material)++; // increment material counter
         num_atoms++; // increment total counter
      }

   }

   // for parallel reduce total atom counter and counts for each material
   #ifdef MPICF

      // variables to store total number of atoms and material associations
      uint64_t total_num_atoms = 0;
      std::vector<uint64_t> material_sum(mp::num_materials,0);

      // now calculate total number of atoms across the system on root process
      MPI_Reduce(&num_atoms,           &total_num_atoms,  1,                 MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&material_numbers[0], &material_sum[0],  mp::num_materials, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);

      // save total num atoms into local variable on root
      num_atoms = total_num_atoms;
      material_numbers = material_sum;

   #endif

   // calculate inverse number of atoms
   const double inverse_num_atoms = 100.0 / double(num_atoms);

	// Output composition to log file
	for(int mat=0;mat<mp::num_materials;mat++){
      zlog << zTs() << "Material " << mat+1 << " " << mp::material[mat].name <<
                       " makes up " << double(material_numbers[mat]) * inverse_num_atoms <<
                       " % of all atoms ( " << material_numbers[mat] <<
                       " atoms )" << std::endl;
   }

	return;

}

} // end of internal namespace

} // end of create namespace
