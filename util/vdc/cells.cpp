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

      if(vdc::verbose) std::cout << "Initialising cell data and allocating atoms to cells..." << std::flush;

      double atoms_min[3] = { 1.0e123, 1.0e123, 1.0e123 };
      double atoms_max[3] = { 0.0, 0.0, 0.0 };

      const double cell_size_x = vdc::cell_size[0]; // cell size (angstroms)
      const double cell_size_y = vdc::cell_size[1]; // cell size (angstroms)
      const double cell_size_z = vdc::cell_size[2]; // cell size (angstroms)

      // calculate minima and maxima
      for(size_t i=0; i < vdc::sliced_atoms_list.size(); i++){

         // get atom ID
         unsigned int atom = vdc::sliced_atoms_list[i];

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
      unsigned int nx = ceil( (atoms_max[0] - atoms_min[0])/cell_size_x );
      unsigned int ny = ceil( (atoms_max[1] - atoms_min[1])/cell_size_y );
      unsigned int nz = ceil( (atoms_max[2] - atoms_min[2])/cell_size_z );

      // check for zero cell size and ensure a minimum of 1 cell in x,y,z
      if( nx == 0 ) nx = 1;
      if( ny == 0 ) ny = 1;
      if( nz == 0 ) nz = 1;

      // save in vdc namespace for calculating newlines for gnuplot compatible 3d data
      vdc::nx_cells = nx;
      vdc::ny_cells = ny;
      vdc::nz_cells = nz;

      unsigned int num_cells[3] = { nx, ny, nz };

      // allocate storage for cell coordinates and cell magnetization
      vdc::total_cells = num_cells[0] * num_cells[1] * num_cells[2];

      // total number of materials + 1
      const unsigned int tmid = 1+vdc::materials.size();

      // arrays to store initial cell corrdinates
      std::vector<double> init_cell_coords(3*total_cells, 0.0);

      // allocate cell memory in 3D and store 1D cell id
      int cell = 0;
      std::vector<std::vector<std::vector<int> > > supercell_array;

      supercell_array.resize(num_cells[0]);
      for(unsigned int i=0;i<num_cells[0];++i){
         supercell_array[i].resize(num_cells[1]);
         for(unsigned int j=0;j<num_cells[1];++j){
            supercell_array[i][j].resize(num_cells[2]);
            // store cell coordinates
            for(unsigned int k=0; k<num_cells[2]; ++k){

               // associate cell with position i,j,k
               supercell_array[i][j][k]=cell;

               // calculate cell coordinates
               init_cell_coords[3*cell + 0] = (double(i) + 0.0) * cell_size_x;
               init_cell_coords[3*cell + 1] = (double(j) + 0.0) * cell_size_y;
               init_cell_coords[3*cell + 2] = (double(k) + 0.0) * cell_size_z;

               // increment cell number
               cell++;
            }
         }
      }

      // Allocate storage for cell id for each atom
      vdc::atom_cell_id.resize(vdc::sliced_atoms_list.size(),0);

      // Determine number of cells in x,y,z (unused)
      // const unsigned int d[3] = { num_cells[0], num_cells[1], num_cells[2] };

      // calculate number of atoms in each cell
      std::vector<unsigned int> init_num_atoms_in_cell(vdc::total_cells);

      //------------------------------------------------------------------------
      // Assign atoms to cells
      //------------------------------------------------------------------------
      for(size_t i=0; i < vdc::sliced_atoms_list.size(); i++){

         // get atom ID
         unsigned int atom = vdc::sliced_atoms_list[i];

         // temporary for atom coordinates
         double c[3] = { vdc::coordinates[3*atom+0] - atoms_min[0],
                         vdc::coordinates[3*atom+1] - atoms_min[1],
                         vdc::coordinates[3*atom+2] - atoms_min[2] };

         // Determine supercell coordinates for atom (rounding down)
         int scc[3] = { int( c[0] / cell_size_x ),
                        int( c[1] / cell_size_y ),
                        int( c[2] / cell_size_z ) };

         // check super cell coordinates are in range
         /*bool error = false;
         for(int j=0; j<3; j++) if(scc[j] >= d[j]) error=true;
         if(error){
            std::cout << "Error! super cell coordinates " << scc[0] << "," << scc[1] << "," << scc[2] << " are out of maximum compute range " << d[0] << "," << d[1] << "," << d[2] << std::endl;
            exit(1);
         }*/

         // Assign atom to cell
         int cellid = supercell_array[scc[0]][scc[1]][scc[2]];
         vdc::atom_cell_id[atom] = cellid;

         // accumulate number of atoms in each cell
         init_num_atoms_in_cell[cellid]++;

      }

      // accumulate total number of cells with atoms
      unsigned num_cells_with_atoms = 0;

      // determine old and new cell numbers (default is number of cells, out of range)
      std::vector<unsigned int> new_cell_number(total_cells, total_cells);

      for( unsigned int cell = 0; cell < total_cells; cell++){
         if( init_num_atoms_in_cell[cell] > 0){

            // save new cell number for this cell
            new_cell_number[cell] = num_cells_with_atoms;

            // save cell coordinates
            vdc::cell_coords.push_back( init_cell_coords[3*cell + 0] );
            vdc::cell_coords.push_back( init_cell_coords[3*cell + 1] );
            vdc::cell_coords.push_back( init_cell_coords[3*cell + 2] );

            // save number of atoms in final cell array
            vdc::num_atoms_in_cell.push_back(init_num_atoms_in_cell[cell]);

            // increment number of cells;
            num_cells_with_atoms++;

         }
      }

      // determine new cell numbers for each atom
      for(size_t i=0; i < vdc::sliced_atoms_list.size(); i++){
         const unsigned int old_cell = vdc::atom_cell_id[i];
         const unsigned int new_cell = new_cell_number[old_cell];

         // check that new cell number is less than total number of cells with atoms
         if(new_cell < num_cells_with_atoms){
            vdc::atom_cell_id[i] = new_cell;
         }
         else{
            std::cerr << "Error - new cell ID " << new_cell << " is larger than number of cells " << num_cells_with_atoms << " with atoms in cell initialisation." << std::endl;
         }

      }

      // Now redfine total number of atoms
      vdc::total_cells = num_cells_with_atoms;

      vdc::cell_magnetization.resize(vdc::total_cells);
      for(unsigned int i = 0; i < vdc::total_cells; i++){
         vdc::cell_magnetization[i].resize(tmid);
         // stored as mx, my, mz, |m| sets
         for(unsigned int j=0; j< tmid; j++) vdc::cell_magnetization[i][j].resize(4,0.0); // initialise cell magnetization to zero
      }

      if(vdc::verbose) std::cout << " done!" << std::endl;

      return;

   }

   void output_cell_file(unsigned int spin_file_id){

      // output informative message to user
      if(vdc::verbose) std::cout << "   Calculating cell magnetization for " << vdc::cell_magnetization.size() << " cells" << std::endl;

      // total number of materials + 1
      const unsigned int tmid = 1+vdc::materials.size();

      // initialise magnetization to zero
      for(unsigned int cell = 0; cell < vdc::total_cells; cell++){
         for(unsigned int m = 0; m < tmid; m++){
            for(int e = 0; e < 4; e++){
               vdc::cell_magnetization[cell][m][e] = 0.0;
            }
         }
      }

      // calculate cell magnetizations in 1D
      for(size_t i=0; i < vdc::sliced_atoms_list.size(); i++){

         // get atom ID
         unsigned int atom = vdc::sliced_atoms_list[i];

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
      for(unsigned int cell = 0; cell < vdc::total_cells; cell++){
         for(unsigned int m = 0; m < tmid; m++){
            const double mx = vdc::cell_magnetization[cell][m][0];
            const double my = vdc::cell_magnetization[cell][m][1];
            const double mz = vdc::cell_magnetization[cell][m][2];
            const double mm = vdc::cell_magnetization[cell][m][3];

            const double norm = sqrt(mx*mx + my*my + mz*mz);

            // calculate inverse norm if norm is greater than 1e-9, otherwise zero
            const double inorm = norm < 1.0e-9 ? 0.0 : 1.0/norm;

            vdc::cell_magnetization[cell][m][0] = mx*inorm;
            vdc::cell_magnetization[cell][m][1] = my*inorm;
            vdc::cell_magnetization[cell][m][2] = mz*inorm;

            // set magnetization of final cell to actual magnetization in mu_B
            if(m == tmid -1) vdc::cell_magnetization[cell][m][3] = norm; // mu_B
            // Otherwise normalise for material magnetization
            else mm < 1.0e-9 ? 0.0 : vdc::cell_magnetization[cell][m][3] = norm/mm; // m/m_s

         }
      }

      // only output cells file if needed
      if(vdc::cellsf){

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

      // write header
      ofile << "#cellx\tcelly\tcellz\t";
      for( unsigned int m = 0; m < tmid-1; m++) ofile << "m"<< m << "mx\t" << "m"<< m << "my\t" << "m"<< m << "mz\t" <<  "m"<< m << "|m|\t";
      ofile << "tot_mx\t" << "tot_my\t" << "tot_mz\t" << "tot_m(muB)\t";
      ofile << "num_atoms" << std::endl;

      for( unsigned int cell = 0; cell < total_cells; cell++){
         ofile << vdc::cell_coords[3*cell + 0] << "\t" << vdc::cell_coords[3*cell + 1] << "\t" << vdc::cell_coords[3*cell + 2] << "\t";
         for( unsigned int m = 0; m < tmid; m++){
            ofile << vdc::cell_magnetization[cell][m][0] << "\t" << vdc::cell_magnetization[cell][m][1] << "\t" << vdc::cell_magnetization[cell][m][2] << "\t" << vdc::cell_magnetization[cell][m][3] << "\t";
         }
         ofile << num_atoms_in_cell[cell];
         ofile << "\n";

         // output new lines after each row of x,y,z for gnuplot compatible data (to be fixed)
         //bool nlx = ( (cell+1) % vdc::nx_cell == 0);
         //bool nlx = ( (cell+1) % vdc::ny_cell == 0);
         //bool nlx = ( (cell+1) % vdc::nz_cell == 0);
         //if((cell+1)%1000 == 0 ) ofile << "\n"; // gnuplot format
      }

      ofile.close();

      std::cout << "done!" << std::endl;

      }

      return;

   }

} // end of namespace vdc
