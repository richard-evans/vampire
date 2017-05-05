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

   //---------------------------------------------------------------------------
   // Function to collate atomistic data into cells to reduce the data density
   //---------------------------------------------------------------------------
   void initialise_cells(){

      double atoms_min[3] = { 1.0e123, 1.0e123, 1.0e123 };
      double atoms_max[3] = { 0.0, 0.0, 0.0 };

      const double cell_size = 10.0; // cell size (angstroms)

      // calculate minima and maxima
      for(unsigned int atom = 0; atom < vdc::num_atoms; atom++){

         const double x = vdc::coordinates[3*atom+0];
         const double y = vdc::coordinates[3*atom+1];
         const double z = vdc::coordinates[3*atom+2];

         if(x < atoms_min[0]) atoms_min[0] = x;
         if(y < atoms_min[1]) atoms_min[1] = y;
         if(z < atoms_min[2]) atoms_min[2] = z;

         if(x > atoms_max[0]) atoms_max[0] = x;
         if(y > atoms_max[1]) atoms_max[1] = y;
         if(z > atoms_max[2]) atoms_max[2] = z;

      }

      // calculate number of cells
      unsigned int nx = ceil( (atoms_max[0] - atoms_min[0])/cell_size );
      unsigned int ny = ceil( (atoms_max[1] - atoms_min[1])/cell_size );
      unsigned int nz = ceil( (atoms_max[2] - atoms_min[2])/cell_size );

      // save in vdc namespace for calculating newlines for gnuplot compatible 3d data
      vdc::nx_cells = nx;
      vdc::ny_cells = ny;
      vdc::nz_cells = nz;

      unsigned int num_cells[3] = { nx, ny, nz };

      // allocate storage for cell coordinates and cell magnetization
      vdc::total_cells = num_cells[0] * num_cells[1] * num_cells[2];

      // total number of materials + 1
      const unsigned int tmid = 1+vdc::materials.size();

      vdc::cell_coords.resize(3*total_cells, 0.0);
      vdc::cell_magnetization.resize(total_cells);
      for(int i = 0; i < vdc::total_cells; i++){
         vdc::cell_magnetization[i].resize(tmid);
         // stored as mx, my, mz, |m| sets
         for(int j=0; j< tmid; j++) vdc::cell_magnetization[i][j].resize(4,0.0); // initialise cell magnetization to zero
      }

      // allocate cell memory in 3D and store 1D cell id
      int cell = 0;
      std::vector<std::vector<std::vector<int> > > supercell_array;

      supercell_array.resize(num_cells[0]);
      for(int i=0;i<num_cells[0];++i){
         supercell_array[i].resize(num_cells[1]);
         for(int j=0;j<num_cells[1];++j){
            supercell_array[i][j].resize(num_cells[2]);
            // store cell coordinates
            for(int k=0; k<num_cells[2]; ++k){

               // associate cell with position i,j,k
               supercell_array[i][j][k]=cell;

               // calculate cell coordinates
               cell_coords[3*cell + 0] = (double(i) + 0.0) * cell_size;
               cell_coords[3*cell + 1] = (double(j) + 0.0) * cell_size;
               cell_coords[3*cell + 2] = (double(k) + 0.0) * cell_size;

               // increment cell number
               cell++;
            }
         }
      }

      // Allocate storage for cell id for each atom
      vdc::atom_cell_id.resize(vdc::num_atoms,0);

      // Determine number of cells in x,y,z
      const unsigned int d[3] = { num_cells[0], num_cells[1], num_cells[2] };

      // Assign atoms to cells
      for(int atom = 0; atom < vdc::num_atoms; atom++ ){

         // temporary for atom coordinates
         double c[3] = { vdc::coordinates[3*atom+0] - atoms_min[0],
                         vdc::coordinates[3*atom+1] - atoms_min[1],
                         vdc::coordinates[3*atom+2] - atoms_min[2] };

         // Determine supercell coordinates for atom (rounding down)
         int scc[3] = { int( c[0] / cell_size ),
                        int( c[1] / cell_size ),
                        int( c[2] / cell_size ) };

         // Assign atom to cell
         vdc::atom_cell_id[atom] = supercell_array[scc[0]][scc[1]][scc[2]];


      }

      return;

   }

   void output_cell_file(unsigned int spin_file_id){

      // output informative message to user
      if(vdc::verbose) std::cout << "   Calculating cell magnetization for " << vdc::cell_magnetization.size() << " cells" << std::endl;

      // total number of materials + 1
      const unsigned int tmid = 1+vdc::materials.size();

      // initialise magnetization to zero
      for(int cell = 0; cell < vdc::total_cells; cell++){
         for(int m = 0; m < tmid; m++){
            for(int e = 0; e < 4; e++){
               vdc::cell_magnetization[cell][m][e] = 0.0;
            }
         }
      }

      // calculate cell magnetizations in 1D
      for(unsigned int atom = 0; atom < vdc::num_atoms; atom++){

         const unsigned int mat = vdc::type[atom];
         const double mu = vdc::materials[mat].moment;

         const double sx = vdc::spins[3*atom+0];
         const double sy = vdc::spins[3*atom+1];
         const double sz = vdc::spins[3*atom+2];

         const unsigned int cell_id = atom_cell_id[atom];

         vdc::cell_magnetization[cell_id][mat][0] += sx*mu;
         vdc::cell_magnetization[cell_id][mat][1] += sy*mu;
         vdc::cell_magnetization[cell_id][mat][2] += sz*mu;
         vdc::cell_magnetization[cell_id][mat][3] += mu;

         // total magnetization in last set
         vdc::cell_magnetization[cell_id][tmid-1][0] += sx*mu;
         vdc::cell_magnetization[cell_id][tmid-1][1] += sy*mu;
         vdc::cell_magnetization[cell_id][tmid-1][2] += sz*mu;
         vdc::cell_magnetization[cell_id][tmid-1][3] += mu;

      }

      // normalise magnetizations
      for(int cell = 0; cell < vdc::total_cells; cell++){
         for(int m = 0; m < tmid; m++){
            const double mx = vdc::cell_magnetization[cell][m][0];
            const double my = vdc::cell_magnetization[cell][m][1];
            const double mz = vdc::cell_magnetization[cell][m][2];
            const double mm = vdc::cell_magnetization[cell][m][3];

            const double norm = sqrt(mx*mx + my*my + mz*mz);
            const double inorm = 1.0/norm;

            vdc::cell_magnetization[cell][m][0] = mx*inorm;
            vdc::cell_magnetization[cell][m][1] = my*inorm;
            vdc::cell_magnetization[cell][m][2] = mz*inorm;

            // set magnetization of final cell to actual magnetization in mu_B
            if(m == tmid -1) vdc::cell_magnetization[cell][m][3] = norm; // mu_B
            // Otherwise normalise for material magnetization
            else vdc::cell_magnetization[cell][m][3] = norm/mm; // m/m_s

         }
	 }

      // output cells to disk
      std::ofstream ofile;

      // Determine cell file name
      std::stringstream cell_file_sstr;
      cell_file_sstr << "cells-";
      cell_file_sstr << std::setfill('0') << std::setw(8) << spin_file_id;
      cell_file_sstr << ".txt";
      std::string cell_file_name = cell_file_sstr.str();

      // output informative message to user
      std::cout << "   Writing cell file " << cell_file_name << "..." << std::flush;

      ofile.open(cell_file_name.c_str());

      for(int cell = 0; cell < total_cells; cell++){
         ofile << vdc::cell_coords[3*cell + 0] << "\t" << vdc::cell_coords[3*cell + 1] << "\t" << vdc::cell_coords[3*cell + 2] << "\t";
         for(int m = 0; m < tmid; m++){
            ofile << vdc::cell_magnetization[cell][m][0] << "\t" << vdc::cell_magnetization[cell][m][1] << "\t" << vdc::cell_magnetization[cell][m][2] << "\t" << vdc::cell_magnetization[cell][m][3] << "\t";
         }
         ofile << "\n";

         // output new lines after each row of x,y,z for gnuplot compatible data (to be fixed)
         //bool nlx = ( (cell+1) % vdc::nx_cell == 0);
         //bool nlx = ( (cell+1) % vdc::ny_cell == 0);
         //bool nlx = ( (cell+1) % vdc::nz_cell == 0);
         //if((cell+1)%1000 == 0 ) ofile << "\n"; // gnuplot format
      }

      ofile.close();

      std::cout << "done!" << std::endl;

      return;

   }

} // end of namespace vdc
