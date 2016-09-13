//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo 2016. All rights reserved.
//
//   Email: am1808@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

// Vampire headers
#include "cells.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

#include "atoms.hpp"

// cells module headers
#include "internal.hpp"

struct uvec{
  int i;
  int j;
  int k;
};

namespace cells{

   //----------------------------------------------------------------------------
   // Function to initialize cells module
   //----------------------------------------------------------------------------
   void initialize(const double system_dimensions_x,
                   const double system_dimensions_y,
                   const double system_dimensions_z,
                   const double unit_cell_size_x,
                   const double unit_cell_size_y,
                   const double unit_cell_size_z,
                   const std::vector<double>& atom_coords_x,
                   const std::vector<double>& atom_coords_y,
                   const std::vector<double>& atom_coords_z,
                   const std::vector<int>& atom_type_array,
                   const std::vector<int>& atom_cell_array,
                   const int num_atoms
   ){

       //-------------------------------------------------------------------------------------
       // Check for cells calculation enabled, if not do nothing
       //-------------------------------------------------------------------------------------
/*       if(!cells::internal::enabled) return;

       // output informative message
       zlog << zTs() << "Initialising data structures for macro-cell calculation." << std::endl;

       // check for prior initialisation
       if(cells::internal::initialised){
          zlog << zTs() << "Warning: Localised temperature calculation already initialised. Continuing." << std::endl;
          return;
       }
*/
      // check calling of routine if error checking is activated
      if(err::check==true) std::cout << "cells::initialise has been called" << std::endl;

      //-------------------------------------------------------------------------------------
      // Define variable needed for mag() function
      //-------------------------------------------------------------------------------------

      cells::internal::num_atoms       = num_atoms;
      cells::internal::atom_type_array = atom_type_array;
      cells::atom_cell_array           = atom_cell_array;

      //-------------------------------------------------------------------------------------
      // Calculate number of microcells
      //-------------------------------------------------------------------------------------
      // determine number of stacks in x and y (global)
      unsigned int dx =  static_cast<unsigned int>(ceil((system_dimensions_x+0.01)/cells::macro_cell_size));
      unsigned int dy =  static_cast<unsigned int>(ceil((system_dimensions_y+0.01)/cells::macro_cell_size));
      unsigned int dz =  static_cast<unsigned int>(ceil((system_dimensions_z+0.01)/cells::macro_cell_size));

      cells::num_cells = dx*dy*dz;
      cells::internal::cell_position_array.resize(3*cells::num_cells);

      //std::cout << " variable cells::num_cells = " << cells::num_cells << std::endl;
      zlog << zTs() << "Macrocell size = " << cells::macro_cell_size << " Angstroms" << std::endl;
      zlog << zTs() << "Macrocells in x,y,z: " << dx << "\t" << dy << "\t" << dz << std::endl;
      zlog << zTs() << "Total number of macrocells: " << cells::num_cells << std::endl;
      zlog << zTs() << "Memory required for macrocell arrays: " << 80.0*double(cells::num_cells)/1.0e6 << " MB" << std::endl;

      //---------------------------------------------------
      // Determine which atoms belong to which cell
      //---------------------------------------------------

      int ncx = dx; // temporary variables for readability
      int ncy = dy;
      int ncz = dz;

      // Set cell and stack counters
      int cell=0;
      //int stack=0;

      // allocate temporary array for neighbour list calculation
      std::vector<uvec> cell_list;
      cell_list.reserve(dx*dy*dz);

      // Allocate space for 3D supercell array (ST coordinate system)
      std::vector<std::vector<std::vector<int> > > supercell_array;
      supercell_array.resize(ncx);
      for(int i=0;i<ncx;++i){
         supercell_array[i].resize(ncy);
         for(int j=0;j<ncy;++j){
            supercell_array[i][j].resize(ncz);
            // store cell coordinates
            for(int k=0; k<ncz; ++k){
               // associate cell with position i,j,k
               supercell_array[i][j][k]=cell;
               // increment cell number
               cell++;
            }
         }
      }

      // Determine number of cells in x,y,z
      const int d[3]={ncx,ncy,ncz};
      const double cs[3] = {cells::macro_cell_size, cells::macro_cell_size, cells::macro_cell_size}; // cell size

     // For MPI version, only add local atoms
      #ifdef MPICF
         int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
      #else
         int num_local_atoms = num_atoms;
      #endif

      // Assign atoms to cells
      for(int atom=0;atom<num_local_atoms;atom++){
         // temporary for atom coordinates
         double c[3];
         // convert atom coordinates to st reference frame
         c[0]=atom_coords_x[atom]+0.0001;
         c[1]=atom_coords_y[atom]+0.0001;
         c[2]=atom_coords_z[atom]+0.0001;
         int scc[3]={0,0,0}; // super cell coordinates
         // Determine supercell coordinates for atom (rounding down)
         scc[0]=int(c[0]/cs[0]);
         scc[1]=int(c[1]/cs[1]);
         scc[2]=int(c[2]/cs[2]);
         for(int i=0;i<3;i++){
            // Always check cell in range
            if(scc[i]<0 || scc[i]>= d[i]){
               terminaltextcolor(RED);
               std::cerr << "Error - atom out of supercell range in cell calculation!" << std::endl;
               terminaltextcolor(WHITE);
               #ifdef MPICF
               terminaltextcolor(RED);
               std::cerr << "\tCPU Rank: " << vmpi::my_rank << std::endl;
               terminaltextcolor(WHITE);
               #endif
               terminaltextcolor(RED);
               std::cerr << "\tAtom number:      " << atom << std::endl;
               std::cerr << "\tCell size:        " << cs[0] << "\t" << cs[1] << "\t" << cs[2] << "\t" << std::endl;
               std::cerr << "\tAtom coordinates: " << c[0] << "\t" << c[1] << "\t" << c[2] << "\t" << std::endl;
               std::cerr << "\tReal coordinates: " << atom_coords_x[atom] << "\t" << atom_coords_y[atom] << "\t" << atom_coords_z[atom] << "\t" << std::endl;
               std::cerr << "\tCell coordinates: " << scc[0] << "\t" << scc[1] << "\t" << scc[2] << "\t" << std::endl;
               std::cerr << "\tCell maxima:      " << d[0] << "\t" << d[1] << "\t" << d[2] << std::endl;
               terminaltextcolor(WHITE);
               err::vexit();
            }
         }
         // If no error for range then assign atom to cell
         cells::atom_cell_array[atom] = supercell_array[scc[0]][scc[1]][scc[2]];

      }

      //-------------------------------------------------------------------------------------
      // Determine number of microcells computed locally
      //-------------------------------------------------------------------------------------

      // Resize new cell arrays
      cells::cell_coords_array_x.resize(cells::num_cells,0.0);
      cells::cell_coords_array_y.resize(cells::num_cells,0.0);
      cells::cell_coords_array_z.resize(cells::num_cells,0.0);

      cells::mag_array_x.resize(cells::num_cells,0.0);
      cells::mag_array_y.resize(cells::num_cells,0.0);
      cells::mag_array_z.resize(cells::num_cells,0.0);

      cells::field_array_x.resize(cells::num_cells,0.0);
      cells::field_array_y.resize(cells::num_cells,0.0);
      cells::field_array_z.resize(cells::num_cells,0.0);

      cells::num_atoms_in_cell.resize(cells::num_cells,0);
      cells::volume_array.resize(cells::num_cells,0.0);

      cells::internal::total_moment_array.resize(cells::num_cells,0.0);

      // Now add atoms to each cell as magnetic 'centre of mass'
      for(int atom=0;atom<num_local_atoms;atom++){
          int local_cell=cells::atom_cell_array[atom];
          int type = cells::internal::atom_type_array[atom];
          const double mus = mp::material[type].mu_s_SI;
          // Consider only magnetic elements
          if(mus/(9.274e-24) > 0.5){
              cells::cell_coords_array_x[local_cell]+=atom_coords_x[atom]*mus;
              cells::cell_coords_array_y[local_cell]+=atom_coords_y[atom]*mus;
              cells::cell_coords_array_z[local_cell]+=atom_coords_z[atom]*mus;

              cells::internal::total_moment_array[local_cell]+=mus;
              cells::num_atoms_in_cell[local_cell]++;
          }
      }

      // For MPI sum coordinates from all CPUs
      #ifdef MPICF
          MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&cells::num_atoms_in_cell[0],  cells::num_cells,MPI_INT,MPI_SUM);
          MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&cells::cell_coords_array_x[0],cells::num_cells,MPI_DOUBLE,MPI_SUM);
          MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&cells::cell_coords_array_y[0],cells::num_cells,MPI_DOUBLE,MPI_SUM);
          MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&cells::cell_coords_array_z[0],cells::num_cells,MPI_DOUBLE,MPI_SUM);
          MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&cells::internal::total_moment_array[0],cells::num_cells,MPI_DOUBLE,MPI_SUM);
      #endif

      // Used to calculate magnetisation in each cell. Poor approximation when unit cell size ~ system size.
      const double atomic_volume = unit_cell_size_x*unit_cell_size_y*unit_cell_size_z/cells::num_atoms_in_unit_cell;

     // Now find mean coordinates via magnetic 'centre of mass'
     for(int local_cell=0;local_cell<cells::num_cells;local_cell++){
         if(cells::num_atoms_in_cell[local_cell]>0){
            //std::cout << cells::num_atoms_in_cell[local_cell] << " in " << local_cell << std::endl;
            cells::cell_coords_array_x[local_cell] = cells::cell_coords_array_x[local_cell]/(cells::internal::total_moment_array[local_cell]);
            cells::cell_coords_array_y[local_cell] = cells::cell_coords_array_y[local_cell]/(cells::internal::total_moment_array[local_cell]);
            cells::cell_coords_array_z[local_cell] = cells::cell_coords_array_z[local_cell]/(cells::internal::total_moment_array[local_cell]);
            cells::volume_array[local_cell] = double(cells::num_atoms_in_cell[local_cell])*atomic_volume;
         }
     }

     cells::atom_in_cell_coords_array_x.resize(cells::num_cells,std::vector<double>(num_local_atoms));
     cells::atom_in_cell_coords_array_y.resize(cells::num_cells,std::vector<double>(num_local_atoms));
     cells::atom_in_cell_coords_array_z.resize(cells::num_cells,std::vector<double>(num_local_atoms));

     cells::index_atoms_array.resize(cells::num_cells,std::vector<int>(num_local_atoms));

      //Set number of atoms in cell to zero
      for(int cell=0;cell<cells::num_cells;cell++){
         cells::num_atoms_in_cell[cell]=0;
      }

      for(int i=0;i<cells::num_cells;i++){
         for(int j=0;j<num_local_atoms;j++){
            cells::atom_in_cell_coords_array_x[i][j]=0.0;
            cells::atom_in_cell_coords_array_y[i][j]=0.0;
            cells::atom_in_cell_coords_array_z[i][j]=0.0;
            cells::index_atoms_array[i][j]=0;
         }
      }

      // Now re-update num_atoms in cell for local atoms only
      for(int atom=0;atom<num_local_atoms;atom++){
         int local_cell=cells::atom_cell_array[atom];
         int type = cells::internal::atom_type_array[atom];
         const double mus = mp::material[type].mu_s_SI;
         // Consider only magnetic elements
         if(mus/(9.274e-24) > 0.5){
            cells::atom_in_cell_coords_array_x[local_cell][cells::num_atoms_in_cell[local_cell]]=atom_coords_x[atom];
            cells::atom_in_cell_coords_array_y[local_cell][cells::num_atoms_in_cell[local_cell]]=atom_coords_y[atom];
            cells::atom_in_cell_coords_array_z[local_cell][cells::num_atoms_in_cell[local_cell]]=atom_coords_z[atom];
            cells::index_atoms_array[local_cell][cells::num_atoms_in_cell[local_cell]]=atom;
            cells::num_atoms_in_cell[local_cell]++;
         }
      }

      // Calculate number of local cells
      for(int cell=0;cell<cells::num_cells;cell++){
         if(cells::num_atoms_in_cell[cell]!=0){
            //std::cout << cells::num_atoms_in_cell[cell] << " in " << cell << std::endl;
            cells::local_cell_array.push_back(cell);
            cells::num_local_cells++;
            cells::volume_array[cell] = double(cells::num_atoms_in_cell[cell])*atomic_volume;
         }
      }

      zlog << zTs() << "Number of local macrocells on rank " << vmpi::my_rank << ": " << cells::num_local_cells << std::endl;

      // Set initialised flag
      cells::internal::initialised=true;

      // Precalculate cell magnetisation
      cells::mag();

      return;

   }

} // end of cells namespace
