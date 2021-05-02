//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>

// program header
#include "vdc.hpp"

namespace vdc{

//------------------------------------------------------------------------------
// Function to output spins-XXXX.vtu files compatible with paraview
//------------------------------------------------------------------------------
void output_vtk_file(unsigned int spin_file_id){

   // Set VTK file name
	std::stringstream vtk_file_sstr;
	vtk_file_sstr << "spins-";
	vtk_file_sstr << std::setfill('0') << std::setw(8) << spin_file_id;
	vtk_file_sstr << ".vtu";
	std::string vtk_file = vtk_file_sstr.str();

   // open vtk file
   std::ofstream vtkfile;
   vtkfile.open(vtk_file.c_str());

   // output informative message to user
   if(vdc::verbose) std::cout << "   Writing VTK file " << vtk_file << "..." << std::flush;

   const double scx = vdc::system_centre[0];
   const double scy = vdc::system_centre[1];
   const double scz = vdc::system_centre[2];

   vtkfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
   vtkfile << "<UnstructuredGrid>" << std::endl;
   vtkfile << "   <Piece NumberOfPoints=\"" << vdc::sliced_atoms_list.size() << "\" NumberOfCells=\"0\">" << std::endl;
   vtkfile << "      <Points>" << std::endl;
   vtkfile << "         <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;


	for(size_t i=0; i < vdc::sliced_atoms_list.size(); i++){

		// get atom ID
		unsigned int atom = vdc::sliced_atoms_list[i];

      vtkfile << "            " << coordinates[3*atom+0]-scx << " " << coordinates[3*atom+1]-scy << " " << coordinates[3*atom+2]-scz << std::endl;

   }
   vtkfile << "         </DataArray>" << std::endl;
   vtkfile << "      </Points>" << std::endl;
   vtkfile << "      <PointData Vectors=\"vector\">" << std::endl;
   vtkfile << "         <DataArray type=\"Float32\" Name=\"spin\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
	for(size_t i=0; i < vdc::atoms_list.size(); i++){

		// get atom ID
		unsigned int atom = vdc::atoms_list[i];

      vtkfile << spins[3*atom+0] << " " << spins[3*atom+1] << " " << spins[3*atom+2] << " ";
   }
   vtkfile << std::endl;
   vtkfile << "         </DataArray>" << std::endl;
   vtkfile << "         <DataArray type=\"Float32\" Name=\"moment\" format=\"ascii\">" << std::endl;
	for(size_t i=0; i < vdc::sliced_atoms_list.size(); i++){

		// get atom ID
		unsigned int atom = vdc::sliced_atoms_list[i];
      // get atom type
      int type_id = vdc::type[atom];
      const double moment = materials[type_id].moment;
      vtkfile << moment << " ";
   }
   vtkfile << std::endl;
   vtkfile << "         </DataArray>" << std::endl;
   vtkfile << "         <DataArray type=\"Int32\" Name=\"ID\" format=\"ascii\">" << std::endl;
	for(size_t i=0; i < vdc::sliced_atoms_list.size(); i++){

		// get atom ID
		unsigned int atom = vdc::sliced_atoms_list[i];
      // get atom type
      int type_id = vdc::type[atom];
      vtkfile << type_id << " ";
   }
   vtkfile << std::endl;
   vtkfile << "         </DataArray>" << std::endl;
   vtkfile << "      </PointData>" << std::endl;
   vtkfile << "      <Cells>" << std::endl;
   vtkfile << "         <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
   vtkfile << "         </DataArray>" << std::endl;
   vtkfile << "         <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
   vtkfile << "         </DataArray>" << std::endl;
   vtkfile << "         <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
   vtkfile << "         </DataArray>" << std::endl;
   vtkfile << "      </Cells>" << std::endl;
   vtkfile << "   </Piece>" << std::endl;
   vtkfile << "</UnstructuredGrid>" << std::endl;
   vtkfile << "</VTKFile>" << std::endl;

   // output informative message to user
   if(vdc::verbose) std::cout << "done!" << std::endl;

   return;

}

} // end of namespace vdc
