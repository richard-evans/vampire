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

#include "dipole.hpp"

#include "atoms.hpp"

// cells module headers
#include "internal.hpp"

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
                   const std::vector<int>& atom_cell_id_array,
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
      cells::atom_cell_id_array        = atom_cell_id_array;

      //-------------------------------------------------------------------------------------
      // Calculate number of microcells
      //-------------------------------------------------------------------------------------
      // determine number of stacks in x and y (global)
      unsigned int dx =  static_cast<unsigned int>(ceil((system_dimensions_x+0.01)/cells::macro_cell_size));
      unsigned int dy =  static_cast<unsigned int>(ceil((system_dimensions_y+0.01)/cells::macro_cell_size));
      unsigned int dz =  static_cast<unsigned int>(ceil((system_dimensions_z+0.01)/cells::macro_cell_size));

      std::cout << "\n macro cell size = " << double( cells::macro_cell_size ) << "\n" << std::flush;

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
      //std::vector<uvec> cell_list;
      //cell_list.reserve(dx*dy*dz);

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
               std::cout << cell << "\t" << i << "\t" << j << "\t" << k << std::endl;
               std::cout << std::flush;
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
               std::cerr << "Error - atom out of supercell range in dipolar field calculation!" << std::endl;
               terminaltextcolor(WHITE);
               #ifdef MPICF
               terminaltextcolor(RED);
               std::cerr << "\tCPU Rank: " << vmpi::my_rank << std::endl;
               terminaltextcolor(WHITE);
               #endif
               terminaltextcolor(RED);
               std::cerr << "\tAtom number:      " << atom << std::endl;
               std::cerr << "\tAtom coordinates: " << c[0] << "\t" << c[1] << "\t" << c[2] << "\t" << std::endl;
               std::cerr << "\tReal coordinates: " << atom_coords_x[atom] << "\t" << atom_coords_y[atom] << "\t" << atom_coords_z[atom] << "\t" << std::endl;
               std::cerr << "\tCell coordinates: " << scc[0] << "\t" << scc[1] << "\t" << scc[2] << "\t" << std::endl;
               std::cerr << "\tCell maxima:      " << d[0] << "\t" << d[1] << "\t" << d[2] << std::endl;
               terminaltextcolor(WHITE);
               err::vexit();
            }
         }
         // If no error for range then assign atom to cell
         cells::atom_cell_id_array[atom] = supercell_array[scc[0]][scc[1]][scc[2]];

         //std::cout << "atom_cell_id_array " << cells::atom_cell_id_array[atom] << " atom " << atom << "\t" << scc[0] << "\t" << scc[1] << "\t" << scc[2] << std::endl;
         //std::cout << std::flush;
         fprintf(stderr,"atom_cell_id_array %d atom %d\t%d\t%d\t%d on rank %d\n",cells::atom_cell_id_array[atom],atom,scc[0],scc[1],scc[2],vmpi::my_rank);
         std::cerr << std::flush;

      }
      MPI::COMM_WORLD.Barrier();

      //-------------------------------------------------------------------------------------
      // Determine number of microcells computed locally
      //-------------------------------------------------------------------------------------

      // Resize new cell arrays

      cells::pos_and_mom_array.resize(4*cells::num_cells,0.0);

      /*cells::cell_coords_array_x.resize(cells::num_cells,0.0);
      cells::cell_coords_array_y.resize(cells::num_cells,0.0);
      cells::cell_coords_array_z.resize(cells::num_cells,0.0);*/

      cells::mag_array_x.resize(cells::num_cells,0.0);
      cells::mag_array_y.resize(cells::num_cells,0.0);
      cells::mag_array_z.resize(cells::num_cells,0.0);

      cells::field_array_x.resize(cells::num_cells,0.0);
      cells::field_array_y.resize(cells::num_cells,0.0);
      cells::field_array_z.resize(cells::num_cells,0.0);

      cells::num_atoms_in_cell.resize(cells::num_cells,0);
      cells::volume_array.resize(cells::num_cells,0.0);

      cells::internal::total_moment_array.resize(cells::num_cells,0.0);
      //cells::cell_id_array.resize(cells::num_cells);

      // Now add atoms to each cell as magnetic 'centre of mass'
      for(int atom=0;atom<num_local_atoms;atom++){
         int local_cell=cells::atom_cell_id_array[atom];
         //int type = cells::internal::atom_type_array[atom];
         int type = atoms::type_array[atom];
         const double mus = mp::material[type].mu_s_SI;
         // Consider only magnetic elements
         if(mus/(9.274e-24) > 0.5){
            /*cells::cell_coords_array_x[local_cell]+=atom_coords_x[atom]*mus;
            cells::cell_coords_array_y[local_cell]+=atom_coords_y[atom]*mus;
            cells::cell_coords_array_z[local_cell]+=atom_coords_z[atom]*mus;
            cells::internal::total_moment_array[local_cell]+=mus;*/

            cells::pos_and_mom_array[4*local_cell+0] += atom_coords_x[atom]*mus;
            cells::pos_and_mom_array[4*local_cell+1] += atom_coords_y[atom]*mus;
            cells::pos_and_mom_array[4*local_cell+2] += atom_coords_z[atom]*mus;
            cells::pos_and_mom_array[4*local_cell+3] += mus;

            cells::num_atoms_in_cell[local_cell]++;
         }
      }

      // Save local cells index
      for(int local_cell=0;local_cell<cells::num_cells;local_cell++){
         if(cells::num_atoms_in_cell[local_cell]>0){
            // add index of cell only if there are atoms inside
            cells::cell_id_array.push_back(local_cell);
         }
      }

      for(unsigned int local_cell=0; local_cell<cells::cell_id_array.size(); local_cell++){
         //std::cout << "index " << local_cell << " of cell_id_array = " << cells::cell_id_array[local_cell] << std::endl;
         //std::cout << std::flush;
         fprintf(stderr,"index %d of cell_id_array = %d on rank %d\n",local_cell,cells::cell_id_array[local_cell],vmpi::my_rank);
      }
      MPI::COMM_WORLD.Barrier();

      fprintf(stderr,"\n\t--- cells::num_atoms_in_cell.size() = %d cell::num_cells = %d on my_rank = %d\n",cells::num_atoms_in_cell.size(),cells::num_cells,vmpi::my_rank);
      #ifdef MPICF
         MPI_Allreduce(MPI_IN_PLACE, &cells::num_atoms_in_cell[0],     cells::num_atoms_in_cell.size(),    MPI_INT,    MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &cells::pos_and_mom_array[0],     cells::pos_and_mom_array.size(),    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		#endif
      fprintf(stderr,"\n\t+++ cells::num_atoms_in_cell.size() = %d cell::num_cells = %d on my_rank = %d\n",cells::num_atoms_in_cell.size(),cells::num_cells,vmpi::my_rank);

      // Used to calculate magnetisation in each cell. Poor approximation when unit cell size ~ system size.
      const double atomic_volume = unit_cell_size_x*unit_cell_size_y*unit_cell_size_z/cells::num_atoms_in_unit_cell;

      // Now find mean coordinates via magnetic 'centre of mass'
      for(int local_cell=0;local_cell<cells::num_cells;local_cell++){
         if(cells::num_atoms_in_cell[local_cell]>0){
            //std::cout << cells::num_atoms_in_cell[local_cell] << " in " << local_cell << std::endl;
            cells::pos_and_mom_array[4*local_cell+0] = cells::pos_and_mom_array[4*local_cell+0]/(cells::pos_and_mom_array[4*local_cell+3]);
            cells::pos_and_mom_array[4*local_cell+1] = cells::pos_and_mom_array[4*local_cell+1]/(cells::pos_and_mom_array[4*local_cell+3]);
            cells::pos_and_mom_array[4*local_cell+2] = cells::pos_and_mom_array[4*local_cell+2]/(cells::pos_and_mom_array[4*local_cell+3]);

            cells::volume_array[local_cell] = double(cells::num_atoms_in_cell[local_cell])*atomic_volume;
         }
      }

      for(int atom=0; atom<num_local_atoms; atom++){
         /*std::cout << "atom " << atom << " in cell " << cells::atom_cell_id_array[atom] << std::endl;
         std::cout << std::flush;*/
         fprintf(stderr,"atom %d with coordinates %f %f %f in cell %d of coordinates %f %f %f in proc %d\n",atom,atom_coords_x[atom],atom_coords_y[atom],atom_coords_z[atom],cells::atom_cell_id_array[atom],cells::pos_and_mom_array[4*cells::atom_cell_id_array[atom]+0],cells::pos_and_mom_array[4*cells::atom_cell_id_array[atom]+1],cells::pos_and_mom_array[4*cells::atom_cell_id_array[atom]+2],vmpi::my_rank);
         std::cout << std::flush;
      }
      MPI::COMM_WORLD.Barrier();

     cells::atom_in_cell_coords_array_x.resize(cells::num_cells,std::vector<double>(num_local_atoms));
     cells::atom_in_cell_coords_array_y.resize(cells::num_cells,std::vector<double>(num_local_atoms));
     cells::atom_in_cell_coords_array_z.resize(cells::num_cells,std::vector<double>(num_local_atoms));

     cells::index_atoms_array.resize(cells::num_cells,std::vector<int>(num_local_atoms));
     //cells::index_atoms_array1D.reserve(num_local_atoms);
     cells::index_atoms_array1D.resize(num_local_atoms);

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
            //cells::index_atoms_array1D[j]=0;
         }
      }

      // Now re-update num_atoms in cell for local atoms only
      for(int atom=0;atom<num_local_atoms;atom++){
         int local_cell=cells::atom_cell_id_array[atom];
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

      std::cout << "\n\nI'm in here Cells Module\n\n" << std::flush;
      // Reorganise array that associate cell with atoms in 1D structure
      // so to pass it easy in parallel
      int counter_id_at=0;
      for(unsigned int lc=0; lc<cells::cell_id_array.size(); lc++){
         for(int i=0; i<cells::num_atoms_in_cell[cells::cell_id_array[lc]]; i++){
            int atom=cells::index_atoms_array[cells::cell_id_array[lc]][i];
            //cells::index_atoms_array1D.push_back(atom);
            cells::index_atoms_array1D[counter_id_at]=atom;
            //std::cout << "lc " << lc << " cell " << cells::cell_id_array[lc] << " loc_atom " << i << " atom " << atom << "\n" << std::flush;
            fprintf(stderr,"lc %d cell %d with coordinates %f %f %f local_atom %d atom %d\t my_rank %d\n",lc,cells::cell_id_array[lc],cells::pos_and_mom_array[4*lc+0],cells::pos_and_mom_array[4*lc+1],cells::pos_and_mom_array[4*lc+2],i,atom,vmpi::my_rank);
            counter_id_at++;
         }
      }
      std::cout << "cells::cell_id_array.size() = " << cells::cell_id_array.size() << "\n" << std::flush;
      std::cout << "\n--->cells::index_atoms_array1D.size() = " << cells::index_atoms_array1D.size() << "\n\n" << std::flush;

      MPI::COMM_WORLD.Barrier();

      #ifdef MPICF
         for(int cpu=0;cpu<vmpi::num_processors;cpu++){
            if(vmpi::my_rank==cpu){
               for(int i=0;i<cells::num_cells;i++){
                  if(cells::num_atoms_in_cell[i]!=0){
                     for(int j=0;j<num_local_atoms;j++){
                      //std::cout << "my_rank = " << vmpi::my_rank << "\t";
                      //std::cout << cells::index_atoms_array[i][j] << "\t";
                      //std::cout << i << "\t";
                      //std::cout << j << "\t";
                      //std::cout << cells::atom_in_cell_coords_array_x[i][j] << "\t";
                      //std::cout << cells::atom_in_cell_coords_array_y[i][j] << "\t";
                      //std::cout << cells::atom_in_cell_coords_array_z[i][j] << std::endl;
                      //std::cout << std::flush;
                        fprintf(stderr,"my_rank = %d\t %d\t %d \t%d\t %f\t %f\t %f\n",vmpi::my_rank,cells::index_atoms_array[i][j],i,j,cells::atom_in_cell_coords_array_x[i][j],cells::atom_in_cell_coords_array_y[i][j],cells::atom_in_cell_coords_array_z[i][j]);
                     }
                  }
               }
            }
         }
         MPI::COMM_WORLD.Barrier();
      #else
         for(int i=0;i<cells::num_cells;i++){
            if(cells::num_atoms_in_cell[i]!=0){
               for(int j=0;j<num_local_atoms;j++){
               	std::cout << cells::index_atoms_array[i][j] << "\t";
               	std::cout << i << "\t";
               	std::cout << j << "\t";
               	std::cout << cells::atom_in_cell_coords_array_x[i][j] << "\t";
               	std::cout << cells::atom_in_cell_coords_array_y[i][j] << "\t";
               	std::cout << cells::atom_in_cell_coords_array_z[i][j] << std::endl;
                  std::cout << std::flush;
               }
            }
         }
      #endif

      // Calculate number of local cells
      for(int cell=0;cell<cells::num_cells;cell++){
         if(cells::num_atoms_in_cell[cell]!=0){
            //std::cout << cells::num_atoms_in_cell[cell] << " in " << cell << std::endl;
            cells::local_cell_array.push_back(cell);
            cells::num_local_cells++;
            cells::volume_array[cell] = double(cells::num_atoms_in_cell[cell])*atomic_volume;
            //std::cout << "cell " << cell << " has " << cells::num_atoms_in_cell[cell] << " atoms inside\n" << std::flush;
            fprintf(stderr,"cell %d of coords %f %f %f has %d atoms inside on cpu=my_rank= %d\n",cell,cells::pos_and_mom_array[4*cell+0],cells::pos_and_mom_array[4*cell+1],cells::pos_and_mom_array[4*cell+2],cells::num_atoms_in_cell[cell],vmpi::my_rank);
            std::cerr << std::flush;
         }
      }
      std::cout << "cell_id_array.size() = " << cell_id_array.size() << "\n" << std::flush;
      std::cout << "local_cell_array.size() = " << local_cell_array.size() << "\n" << std::flush;
      std::cout << "num_local_cells = " << num_local_cells << "\n" << std::flush;
      for(int lc=0; lc<cells::num_local_cells; lc++){
         #ifdef MPICF
            fprintf(stderr,"%d\t%d\t%f\t%f\t%f\t%e\n",cell_id_array[lc],local_cell_array[lc],pos_and_mom_array[4*cell_id_array[lc]+0],pos_and_mom_array[4*cell_id_array[lc]+1],pos_and_mom_array[4*cell_id_array[lc]+2],pos_and_mom_array[4*cell_id_array[lc]+3]);
         #else
            std::cout << cell_id_array[lc] << "\t";
            std::cout << local_cell_array[lc] << "\t";
            std::cout << pos_and_mom_array[4*cell_id_array[lc]+0] << "\t";
            std::cout << pos_and_mom_array[4*cell_id_array[lc]+1] << "\t";
            std::cout << pos_and_mom_array[4*cell_id_array[lc]+2] << "\t";
            std::cout << pos_and_mom_array[4*cell_id_array[lc]+3] << "\n" << std::flush;
            std::cout << std::endl << std::flush;
         #endif
      }
      MPI::COMM_WORLD.Barrier();

   ///// Reduce to all cpus cells info
   ///#ifdef MPICF
   ///    MPI_Allreduce(MPI_IN_PLACE, &cells::num_atoms_in_cell[0],     cells::num_atoms_in_cell.size(),    MPI_INT,    MPI_SUM, MPI_COMM_WORLD);
   ///    MPI_Allreduce(MPI_IN_PLACE, &cells::pos_and_mom_array[0],     cells::pos_and_mom_array.size(),    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   ///    //MPI_Allreduce(MPI_IN_PLACE, &cells::cell_id_array[0],         cells::cell_id_array.size(),        MPI_INT,    MPI_SUM, MPI_COMM_WORLD);
   ///    //MPI_Allreduce(MPI_IN_PLACE, &cells::atom_cell_id_array[0],    cells::atom_cell_id_array.size(),   MPI_INT,    MPI_SUM, MPI_COMM_WORLD);
   ///    //MPI_Allreduce(MPI_IN_PLACE, &cells::index_atoms_array1D[0],   cells::index_atoms_array1D.size(),  MPI_INT,    MPI_SUM, MPI_COMM_WORLD);

   ///   //MPI_Bcast(&cells::atom_cell_id_array[0],cells::atom_cell_id_array.size(),MPI_INT,0,MPI_COMM_WORLD);
   ///   //MPI_Bcast(&cells::index_atoms_array1D[0],index_atoms_array1D.size(),MPI_INT,0,MPI_COMM_WORLD);
   ///   //MPI_Allgather(&cells::atom_cell_id_array[0],vmpi::num_core_atoms+vmpi::num_bdry_atoms,MPI_INT,&cells::atom_cell_id_array[0],cells::atom_cell_id_array.size(),MPI_INT,MPI_COMM_WORLD);
   ///   //MPI_Allgather(&cells::index_atoms_array1D[0],vmpi::num_core_atoms+vmpi::num_bdry_atoms,MPI_INT,&cells::index_atoms_array1D[0],index_atoms_array1D.size(),MPI_INT,MPI_COMM_WORLD);
   ///   /*if(vmpi::my_rank==0){
   ///      for(int cpu=1; cpu<vmpi::num_processors; cpu++){
   ///         std::cout << "\nI am the root --> cpu = " << vmpi::my_rank << std::endl << std::flush;
   ///         MPI_Send(&cells::atom_cell_id_array[0],cells::atom_cell_id_array.size(),MPI_INT,cpu,101,MPI_COMM_WORLD);
   ///         MPI_Send(&cells::index_atoms_array1D[0],index_atoms_array1D.size(),MPI_INT,cpu,101,MPI_COMM_WORLD);
   ///         std::cout << "root is sending data to slave " << cpu << std::endl << std::flush;
   ///      }
   ///   }
   ///   else{ // I am a slave
   ///      std::cout << "\nI am NOT the root --> cpu = " << vmpi::my_rank << std::endl << std::flush;
   ///      MPI_Recv(&cells::atom_cell_id_array[0],cells::atom_cell_id_array.size(),MPI_INT,0,101,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
   ///      MPI_Recv(&cells::index_atoms_array1D[0],index_atoms_array1D.size(),MPI_INT,0,101,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
   ///      std::cout << vmpi::my_rank << " is receiving data from root\n" << std::flush;
   ///   }*/
   ///   std::cout << std::endl << std::flush;
	///#endif
      std::cout << "\n--->cells::index_atoms_array1D.size() = " << cells::index_atoms_array1D.size() << "\n\n" << std::flush;

      std::vector<double> atom_pos_x(num_local_atoms,0.0);
      std::vector<double> atom_pos_y(num_local_atoms,0.0);
      std::vector<double> atom_pos_z(num_local_atoms,0.0);
      for(int atom=0; atom<num_local_atoms; atom++){
         atom_pos_x[atom]=atom_coords_x[atom];
         atom_pos_y[atom]=atom_coords_y[atom];
         atom_pos_z[atom]=atom_coords_z[atom];
      }

      // Call paralelisation function
      //dipole::send_receive_data(num_local_atoms,atom_pos_x,atom_pos_y,atom_pos_z);

//    if(vmpi::my_rank==0){
//       for(int lc=0; lc<cells::num_local_cells; lc++){
//          /*std::cout << "\t";
//          std::cout << "my_rank==0\t";
//          std::cout << cell_id_array[lc] << "\t";
//          std::cout << pos_and_mom_array[4*cell_id_array[lc]+0] << "\t";
//          std::cout << pos_and_mom_array[4*cell_id_array[lc]+1] << "\t";
//          std::cout << pos_and_mom_array[4*cell_id_array[lc]+2] << "\t";
//          std::cout << pos_and_mom_array[4*cell_id_array[lc]+3] << "\n" << std::flush;*/
//          fprintf(stderr,"\tmyrank ==%d\t%d\t%f\t%f\t%f\t%f\n",vmpi::my_rank,cell_id_array[lc],pos_and_mom_array[4*cell_id_array[lc]+0],pos_and_mom_array[4*cell_id_array[lc]+1],pos_and_mom_array[4*cell_id_array[lc]+2],pos_and_mom_array[4*cell_id_array[lc]+3]);
//       }
//    }
//    else{
//       for(int lc=0; lc<cells::num_local_cells; lc++){
//          /*std::cout << "\t";
//          std::cout << "my_rank = " << vmpi::my_rank << "\t";
//          std::cout << cell_id_array[lc] << "\t";
//          std::cout << pos_and_mom_array[4*cell_id_array[lc]+0] << "\t";
//          std::cout << pos_and_mom_array[4*cell_id_array[lc]+1] << "\t";
//          std::cout << pos_and_mom_array[4*cell_id_array[lc]+2] << "\t";
//          std::cout << pos_and_mom_array[4*cell_id_array[lc]+3] << "\n" << std::flush;*/
//          fprintf(stderr,"\tmyrank = %d\t%d\t%f\t%f\t%f\t%f\n",vmpi::my_rank,cell_id_array[lc],pos_and_mom_array[4*cell_id_array[lc]+0],pos_and_mom_array[4*cell_id_array[lc]+1],pos_and_mom_array[4*cell_id_array[lc]+2],pos_and_mom_array[4*cell_id_array[lc]+3]);
//       }
//    }

      for(int lc=0; lc<cells::num_local_cells; lc++){
            fprintf(stderr,"\tmy_rank=%d\tcell_id_array[lc] %d\tx %f\ty %f\tz %f\tmu %e\n",vmpi::my_rank,cell_id_array[lc],pos_and_mom_array[4*cell_id_array[lc]+0],pos_and_mom_array[4*cell_id_array[lc]+1],pos_and_mom_array[4*cell_id_array[lc]+2],pos_and_mom_array[4*cell_id_array[lc]+3]);
      }


//    for(int lc=0; lc<cells::num_local_cells; lc++){
//       if(vmpi::my_rank==0){
//          fprintf(stderr,"\tmyrank ==%d\t%d\t%f\t%f\t%f\t%f\n",vmpi::my_rank,cell_id_array[lc],pos_and_mom_array[4*cell_id_array[lc]+0],pos_and_mom_array[4*cell_id_array[lc]+1],pos_and_mom_array[4*cell_id_array[lc]+2],pos_and_mom_array[4*cell_id_array[lc]+3]);
//       }
//       else{
//          fprintf(stderr,"\tmyrank = %d\t%d\t%f\t%f\t%f\t%f\n",vmpi::my_rank,cell_id_array[lc],pos_and_mom_array[4*cell_id_array[lc]+0],pos_and_mom_array[4*cell_id_array[lc]+1],pos_and_mom_array[4*cell_id_array[lc]+2],pos_and_mom_array[4*cell_id_array[lc]+3]);
//       }
//    }

      //std::cout << "\n\n\n" << std::flush;
      fprintf(stderr,"\n\n");
      MPI::COMM_WORLD.Barrier();

//    for(int lc=0; lc<cells::num_local_cells; lc++){
//       for(int i=0; i<cells::num_atoms_in_cell[cell_id_array[lc]]; i++){
//          if(vmpi::my_rank==0){
//             fprintf(stderr,"\tmyrank ==%d\tindex_atoms_array1D[%d] = %d\tin cell = %d\n",vmpi::my_rank,lc+i,index_atoms_array1D[i+lc],lc);
//          }
//          else{
//             fprintf(stderr,"\tmyrank = %d\tindex_atoms_array1D[%d] = %d\tin cell = %d\n",vmpi::my_rank,lc+i,index_atoms_array1D[i+lc],lc);
//          }
//       }
//    }

      int counter_cell=0;
//    if(vmpi::my_rank==0){
//       for(int lc=0; lc<cells::num_local_cells; lc++){
//          for(int i=0; i<cells::num_atoms_in_cell[cell_id_array[lc]]; i++){
//             /*std::cout << "\t";
//             std::cout << "my_rank==0\t";
//             std::cout << index_atoms_array1D[i+lc] << "\t";
//             std::cout << "in cell = " << lc;
//             std::cout << "\n" << std::flush;*/
//             fprintf(stderr,"\tmyrank ==%d\tindex_atoms_array1D[%d] = %d\tin cell = %d\n",vmpi::my_rank,lc+i,index_atoms_array1D[i+lc],lc);
//          }
//       }
//    }
//    else{
//       for(int lc=0; lc<cells::num_local_cells; lc++){
//          for(int i=0; i<cells::num_atoms_in_cell[cell_id_array[lc]]; i++){
//             /*std::cout << "\t";
//             std::cout << "my_rank = " << vmpi::my_rank << "\t";
//             std::cout << index_atoms_array1D[i+lc] << "\t";
//             std::cout << "in cell = " << lc;
//             std::cout << "\n" << std::flush;*/
//             fprintf(stderr,"\tmyrank = %d\tindex_atoms_array1D[%d] = %d\tin cell = %d\n",vmpi::my_rank,lc+i,index_atoms_array1D[i+lc],lc);
//          }
//       }
//    }
      fprintf(stderr,"\n\n");


      zlog << zTs() << "Number of local macrocells on rank " << vmpi::my_rank << ": " << cells::num_local_cells << std::endl;

      // Set initialised flag
      cells::internal::initialised=true;

      //// Precalculate cell magnetisation
      ////cells::mag();

      return;

   }

} // end of cells namespace
