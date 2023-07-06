//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins 2016. All rights reserved.
//
//   Email: sj681@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "environment.hpp"
#include "random.hpp"
#include "../cells/internal.hpp"
#include "../micromagnetic/internal.hpp"
#include "micromagnetic.hpp"

// environment module headers
#include "internal.hpp"
#include "sim.hpp"
#include "atoms.hpp"
#include "cells.hpp"
#include "vmpi.hpp"

#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>




 int x_cubes[] = {1,1,1,1,-1,-1,-1,-1};
 int y_cubes[] = {1,1,-1,-1,-1,-1,1,1};
 int z_cubes[] = {1,-1,-1,1,-1,1,1,-1};



using namespace std;

namespace env = environment::internal;

void neighbour(int env, int mm, double overlap);
void split_cells(double x_pos, double y_pos, double z_pos, double current_cell_size);

namespace environment{

   //----------------------------------------------------------------------------
   // Function to initialize environment module
   //----------------------------------------------------------------------------
   void initialize(double system_dimensions_x, double system_dimensions_y, double system_dimensions_z){


       std::cout << "Initialising environment module" << std::endl;


       env::env_field_uv[0] = env::env_field_uv[0]*env::env_field;
       env::env_field_uv[1] = env::env_field_uv[1]*env::env_field;
       env::env_field_uv[2] = env::env_field_uv[2]*env::env_field;

       //double mm_size_x, mm_size_y,mm_size_z;
       //double env_size_x, env_size_y,env_size_z;

       env::shield_shape.resize(env::num_shields+1,"cube");
       env::shield_ms.resize(env::num_shields+1, 1e-21);
       env::shield_Tc.resize(env::num_shields+1, 600);
       env::shield_A.resize(env::num_shields+1);
       for (int i = 0; i < env::num_shields+1; ++i)
        env::shield_A[i].resize(env::num_shields+1, -180);
       env::shield_ku.resize(env::num_shields+1, 0);
       env::shield_alpha.resize(env::num_shields+1, 1.0);
       env::shield_gamma.resize(env::num_shields+1,1.0);
       env::shield_max_x.resize(env::num_shields+1,1000);
       env::shield_max_y.resize(env::num_shields+1,1000);
       env::shield_max_z.resize(env::num_shields+1,1000);
       env::shield_min_x.resize(env::num_shields+1,500);
       env::shield_min_y.resize(env::num_shields+1,500);
       env::shield_min_z.resize(env::num_shields+1,500);
       env::shield_max_cell_size.resize(env::num_shields+1,100);
       env::shield_min_cell_size.resize(env::num_shields+1,10);
       env::shield_Hext_x.resize(env::num_shields+1,0);
       env::shield_Hext_y.resize(env::num_shields+1,0);
       env::shield_Hext_z.resize(env::num_shields+1,0);
       env::pos_or_neg.resize(env::num_shields+1,"pos");
       env::H_strength.resize(env::num_shields+1,0);
       env::initial_spin_x.resize(env::num_shields+1,0);
       env::initial_spin_y.resize(env::num_shields+1,0);
       env::initial_spin_z.resize(env::num_shields+1,1);
       env::random_spins.resize(env::num_shields+1,false);


       env::read_in_shield_info();

       for (int shield=0; shield < env::num_shields; shield++){
         env::shield_Hext_x[shield] = env::shield_Hext_x[shield]*env::H_strength[shield];
         env::shield_Hext_y[shield] = env::shield_Hext_y[shield]*env::H_strength[shield];
         env::shield_Hext_z[shield] = env::shield_Hext_z[shield]*env::H_strength[shield];
      //   std::cout << env::shield_Hext_x[shield] << '\t' << env::shield_Hext_y[shield] << '\t' << env::shield_Hext_z[shield] <<std::endl;
       }
    //   std::cout << "shield info read in " <<std::endl;

       double dx,dy,dz;

       //stores the max and min coordinates of all env cells
       std::vector <double > x_max;
       std::vector <double > y_max;
       std::vector <double > z_max;
       std::vector <double > x_min;
       std::vector <double > y_min;
       std::vector <double > z_min;

       std::vector < std::vector < double> >corners;
       corners.resize(8);
       for (int i = 0; i < 8; ++i) corners[i].resize(3,0.0);
       std::ofstream pfile2;
       pfile2.open("env_cell_positions");
       env::o_file.open("env_output");

       int n_cell = 0;
      // std::cout << env::num_shields <<std::endl;


       for (int shield = 0; shield < env::num_shields; shield++){

         std::vector <double> tmp_x;
         std::vector <double> tmp_y;
         std::vector <double> tmp_z;
         std::vector <double> size_x;
         std::vector <double> size_y;
         std::vector <double> size_z;


         dx = env::shield_max_x[shield] - env::shield_min_x[shield];
         dy = env::shield_max_y[shield] - env::shield_min_y[shield];
         dz = env::shield_max_z[shield] - env::shield_min_z[shield];
      //   std::cout << dx << '\t' << dy << '\t' << dz << '\t' << env::cell_size <<std::endl;
         env::num_cells_x =  static_cast<unsigned int>(ceil((dx+0.01)/env::shield_max_cell_size[shield]));
         env::num_cells_y =  static_cast<unsigned int>(ceil((dy+0.01)/env::shield_max_cell_size[shield]));
         env::num_cells_z =  static_cast<unsigned int>(ceil((dz+0.01)/env::shield_max_cell_size[shield]));

        // std::cout <<"A" <<  env::num_cells_x << '\t' << env::num_cells_y << '\t' << env::num_cells_z <<std::endl;
          //                    std::cin.get();
         //calculates the minimum and maximum positions of each cell for calcualtions of which cell the atomistic atoms are in later
   //      std::cout << shield << '\t' <<  tmp_x.size() << "\t" << env::cell_coords_array_x.size() << "\t" << n_cell <<std::endl;
   //      std::cout << env::num_cells_x  << '\t' << env::num_cells_y  << "\t"  << env::num_cells_z << "\t" << env::shield_max_cell_size[shield]<< "\t" << dx << "\t" << dy << "\t" << dz << std::endl;
       for (int x = 0; x < env::num_cells_x; x++){
          for (int y = 0; y < env::num_cells_y; y++){
             for (int z = 0; z < env::num_cells_z; z++){

              // std::cout << "here" << '\t' << x << '\t' << y << '\t' << z << std::endl;
                double current_cell_size = env::shield_max_cell_size[shield];
              //  std::cout << shield << '\t' << current_cell_size << std::endl;
                double pos_x = current_cell_size*x + env::shield_min_x[shield] + current_cell_size/2.0;// + current_cell_size/2.0;
                double pos_y = current_cell_size*y + env::shield_min_y[shield] + current_cell_size/2.0 ;// + current_cell_size/2.0;
                double pos_z = current_cell_size*z + env::shield_min_z[shield] + current_cell_size/2.0 ;// + current_cell_size/2.0;

            //   std::cout << "here" << '\t' << pos_x << '\t' << pos_y << '\t' << pos_z << std::endl;
                size_x.push_back(current_cell_size);
                size_y.push_back(current_cell_size);
                size_z.push_back(current_cell_size);
                tmp_x.push_back(pos_x);
                tmp_y.push_back(pos_y);
                tmp_z.push_back(pos_z);

              }
            }
          }

          std::cout <<  tmp_x.size() <<std::endl;
          //       std::cout << "created initial cells " <<std::endl;
        //  std::cin.get();

          for (size_t cell = 0; cell < tmp_x.size(); cell ++ ){

        //    std::cout << shield << '\t' << cell << "\txyzdx\t" << tmp_x[cell] << '\t' << tmp_y[cell] << '\t' << tmp_z[cell]  << '\t' << size_x[cell] << "\t" <<  std::endl;
          //  std::cout << "create corners" << std::endl;
            corners = env::calculate_corners(tmp_x[cell],tmp_y[cell],tmp_z[cell],size_x[cell],size_y[cell],size_z[cell]);
            int N = 0;
            for (int corner = 0; corner < 8; corner ++){

              int within = env::in_shield(corners[corner][0] , corners[corner][1] , corners[corner][2],shield);
            //  std::cout << "corners" << '\t' << within << "\t" << corners[corner][0]  << '\t' << corners[corner][1]  << '\t' << corners[corner][2]  << std::endl;
              if (within != 0){
                  N++;
              }
            }

              if (N != 8 && size_x[cell] > 40.0){
                //  std::cout << shield << '\t' << "SPLIT" << std::endl;
              //    std::cout << "splitting cell \t:" << cell << '\t' <<  tmp_x[cell] << '\t' << tmp_y[cell] << '\t' << tmp_z[cell]  << '\t' << size_x[cell] << std::endl;
                for (int n = 0; n <8; n++ ){
                  double new_cell_size = size_x[cell]/2.0;
                  size_x.push_back(new_cell_size);
                  size_y.push_back(new_cell_size);
                  size_z.push_back(new_cell_size);
                  tmp_x.push_back(tmp_x[cell] + x_cubes[n]*new_cell_size/2.0);
                  tmp_y.push_back(tmp_y[cell] + y_cubes[n]*new_cell_size/2.0);
                  tmp_z.push_back(tmp_z[cell] + z_cubes[n]*new_cell_size/2.0);

              //    std::cout << new_cell_size << '\t' << tmp_x[cell] + x_cubes[n]*new_cell_size/2.0 << '\t' << tmp_y[cell] + y_cubes[n]*new_cell_size/2.0 << '\t' << tmp_z[cell] + z_cubes[n]*new_cell_size/2.0 << "\t" << tmp_x.size() << '\t' << cell << std::endl;
                //  std::cout << new_cell_size <<std::endl;
                }
            //      std::cout << "cell split" <<std::endl;
                }
                else if (N != 0 ) {

                //  std::cout << shield << '\t' << "SAVED" << std::endl;
                  env::cell_coords_array_x.push_back(tmp_x[cell]);
                  env::cell_coords_array_y.push_back(tmp_y[cell]);
                  env::cell_coords_array_z.push_back(tmp_z[cell]);
                  env::cell_size_x.push_back(size_x[cell]);
                  env::cell_size_y.push_back(size_y[cell]);
                  env::cell_size_z.push_back(size_z[cell]);
                  env::cell_volume.push_back(size_x[cell]*size_y[cell]*size_z[cell]);
                  env::Ms.push_back(env::shield_ms[shield]*env::cell_volume[n_cell]);
                  //std::cout << env::shield_ms[shield] << '\t' << env::cell_volume[n_cell] << std::endl;
                  env::ku.push_back(-env::shield_ku[shield]/env::cell_volume[n_cell]);
                  x_max.push_back(tmp_x[cell] + size_x[cell]/2.0);
                  x_min.push_back(tmp_x[cell] - size_x[cell]/2.0);
                  y_max.push_back(tmp_y[cell] + size_y[cell]/2.0);
                  y_min.push_back(tmp_y[cell] - size_y[cell]/2.0);
                  z_max.push_back(tmp_z[cell] + size_z[cell]/2.0);
                  z_min.push_back(tmp_z[cell] - size_z[cell]/2.0);
                  env::x_mag_array.push_back(0.0);
                  env::y_mag_array.push_back(0.0);
                  env::z_mag_array.push_back(0.0);
                  env::shield_number.push_back(shield);
                  env::env_cell_is_in_atomistic_region.push_back(0);
                  pfile2 << N << '\t' << n_cell << '\t' <<  tmp_x[cell] << '\t' << tmp_y[cell] << '\t' << tmp_z[cell] << "\t" << size_x[cell] << '\t' << size_y[cell]<< '\t' << size_z[cell] <<  std::endl;
                  n_cell ++;

                }
              //  std::cin.get();
              }
            //  std::cout << "END\t" << shield << '\t' <<  tmp_x.size() << "\t" << env::cell_coords_array_x.size() << "\t" << n_cell <<std::endl;
            //  std::cin.get();
            }


       env::num_cells = n_cell;
    //   std::cout << "Number of environment cells: " << env::num_cells <<  "\t cell size: " <<   env::cell_volume[0] << "\t" << std::endl;

       //convert Ms from input to Ms = ms/V and Ku = ku/V






  //     std::cout << "Identifying environment cells which are atomistic" << std::endl;

       //loops over all atomistic cells to determine if the atomsitic simulation lies within an environment cell

       // std::ofstream pfile;
       // pfile.open("m3.txt");
       //
       // for(int lc=0; lc<cells::num_local_cells; lc++){
       //    int cell = cells::cell_id_array[lc];
       //    pfile << cell << "\t" << cells::pos_and_mom_array[4*cell + 0] << "\t" << cells::pos_and_mom_array[4*cell + 1] << "\t" << cells::pos_and_mom_array[4*cell + 2] <<  "\t" <<  cells::mag_array_x[cell] << '\t' << cells::mag_array_y[cell] << '\t' << cells::mag_array_z[cell] <<std::endl;
       //
       // }

      env::list_env_cell_atomistic_cell.resize(cells::num_cells,0);

       for(int lc=0; lc<cells::num_local_cells; lc++){
          int cell = cells::cell_id_array[lc];
       //   std::cout << cell << '\t' <<std::endl;
          //calcaultes the x,y,z position for each atomsitic cell
          double x,y,z;
          if (env::shift[0] >= 0) x = cells::pos_and_mom_array[4*cell + 0] + env::shift[0];
          else x = cells::pos_and_mom_array[4*cell + 0];
          if (env::shift[1] >= 0) y = cells::pos_and_mom_array[4*cell + 1] + env::shift[1];
          else y = cells::pos_and_mom_array[4*cell + 1];
          if (env::shift[2] >= 0) z = cells::pos_and_mom_array[4*cell + 2] + env::shift[2];
          else z = cells::pos_and_mom_array[4*cell + 2];
          bool in_cell = false;

          //loops over all environment cells to determine which cell this atomistic cell lies within
          for (int env_cell = 0; env_cell < env::num_cells; env_cell++){

             //if atom is within environment
         //   std::cout << x << '\t' << y << '\t' << z << '\t' << x_max[env_cell] << '\t' << y_max[env_cell] << '\t' << z_max[env_cell] <<std::endl;
             if (x < x_max[env_cell] && x > x_min[env_cell] && y < y_max[env_cell] && y > y_min[env_cell] && z < z_max[env_cell] && z > z_min[env_cell]){
                //then the magnetisation of the enciroment cell is a sum of all atomistic atoms in that cell
                //adds the cell to the list of atomistic cells
                env::list_env_cell_atomistic_cell[cell] = env_cell;
                //this cell is an atomistic cell.
                env::env_cell_is_in_atomistic_region[env_cell] = 1;
                env::x_mag_array[env_cell] += cells::mag_array_x[cell];
                env::y_mag_array[env_cell] += cells::mag_array_y[cell];
                env::z_mag_array[env_cell] += cells::mag_array_z[cell];
                in_cell = true;
          //      std::cout << env_cell << '\t' << cell << std::endl;

             }
          }
          if (in_cell == false){
            env::num_cells ++;
            int env_cell = env::num_cells - 1;
            env::cell_coords_array_x.push_back(x);
            env::cell_coords_array_y.push_back(y);
            env::cell_coords_array_z.push_back(z);
            env::cell_size_x.push_back(cells::macro_cell_size_x);
            env::cell_size_y.push_back(cells::macro_cell_size_y);
            env::cell_size_z.push_back(cells::macro_cell_size_z);
            env::cell_volume.push_back(cells::macro_cell_size_x*cells::macro_cell_size_y*cells::macro_cell_size_z);
            env::Ms.push_back(1e-21);
            env::ku.push_back(0);
            env::shield_number.push_back(env::num_shields);
            env::list_env_cell_atomistic_cell[cell] = env_cell;
            //this cell is an atomistic cell.
            env::env_cell_is_in_atomistic_region.push_back(1);
            env::x_mag_array.push_back(cells::mag_array_x[cell]);
            env::y_mag_array.push_back(cells::mag_array_y[cell]);
            env::z_mag_array.push_back(cells::mag_array_z[cell]);
            pfile2 << 234 << '\t' << env_cell << '\t' <<  x << '\t' << y << '\t' << z << "\t" << cells::macro_cell_size_x << '\t' << cells::macro_cell_size_y << '\t' << cells::macro_cell_size_z <<  std::endl;
          }
       }

    //   std::cout << "end of identification" << std::endl;
       //resize arrays
       env::bias_field_x.resize(env::num_cells,0.0);
       env::bias_field_y.resize(env::num_cells,0.0);
       env::bias_field_z.resize(env::num_cells,0.0);

       env::dipole_field_x.resize(env::num_cells,0.0);
       env::dipole_field_y.resize(env::num_cells,0.0);
       env::dipole_field_z.resize(env::num_cells,0.0);

       env::neighbour_list_start_index.resize(env::num_cells,0.0);
       env::neighbour_list_end_index.resize(env::num_cells,0.0);

       env::env_cell_is_in_atomistic_region.resize(env::num_cells,0);

       // For MPI version, only add local atoms
        #ifdef MPICF
           int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
        #else
           int num_local_atoms = atoms::num_atoms;
        #endif

       environment_field_x.resize(cells::num_cells,0.0);
       environment_field_y.resize(cells::num_cells,0.0);
       environment_field_z.resize(cells::num_cells,0.0);
       atomistic_environment_field_x.resize(num_local_atoms,0.0);
       atomistic_environment_field_y.resize(num_local_atoms,0.0);
       atomistic_environment_field_z.resize(num_local_atoms,0.0);

    //   std::cin.get();
  //  std::cout << "resized arrays" << std::endl;

        // std::ofstream pfile2;
        // pfile2.open("m2.txt");
       //
       // for (int cell = 0; cell < env::num_cells; cell++){
       //    pfile2 << cell << "\t" << env::cell_coords_array_x[cell] << "\t" << env::cell_coords_array_y[cell] << "\t" << env::cell_coords_array_z[cell] <<  "\t" <<  env::x_mag_array[cell] << '\t' << env::y_mag_array[cell] << '\t' << env::z_mag_array[cell] <<std::endl;
       //  }

       //calcualtes the neighbour lists for each cell.
       //if cells are neighbours add them to the neighbour list array_index
       //each cell has a start and end index index. the neighbours for each cell are within these index

   //    std::cout << "Calculating neighbour list for environment cells" << "\t" << env::num_cells <<env::cell_coords_array_x.size() << "\t" << env::cell_size_x.size() << '\t' <<  env::neighbour_list_end_index.size() << std::endl;
       // This is massively inefficient for large numbers of cells
       // better to store x,y,z cell associations and calculate neighbours directly - RE

       int array_index = 0;
       for (int celli = 0; celli < env::num_cells; celli ++){
          double xi = env::cell_coords_array_x[celli];
          double yi = env::cell_coords_array_y[celli];
          double zi = env::cell_coords_array_z[celli];
          env::neighbour_list_start_index[celli] = array_index;
          int b = 0;
          for (int cellj = 0; cellj < env::num_cells; cellj ++){
             int a = 0;
             double xj = env::cell_coords_array_x[cellj];
             double yj = env::cell_coords_array_y[cellj];
             double zj = env::cell_coords_array_z[cellj];
             double dx = sqrt((xi - xj)*(xi - xj));
             double dy = sqrt((yi - yj)*(yi - yj));
             double dz = sqrt((zi - zj)*(zi - zj));

             //calcualtes how many sides of the cube are touching neaarest neighbours have 1.
             if (dx == env::cell_size_x[celli]/2.0 +  env::cell_size_x[cellj]/2.0 && dy ==0 && dz == 0) a++;
             if (dy == env::cell_size_y[celli]/2.0 +  env::cell_size_y[cellj]/2.0 && dz ==0 && dx == 0) a++;
             if (dz == env::cell_size_z[celli]/2.0 +  env::cell_size_z[cellj]/2.0 && dy ==0 && dx == 0) a++;
             if (a == 1){
                b++;
                env::neighbour_list_array.push_back(cellj);                                        //if the interaction is non zero add the cell to the neighbourlist
                env::neighbour_list_end_index[celli] = array_index;                                //the end index is updated for each cell so is given the value for the last cell.
                array_index ++;
          //      std::cout << celli << '\t' << cellj << std::endl;
             }
          }
       }

    //   std::cout << "End of enighbou list calculation" << "\t" << env::num_cells <<env::cell_coords_array_x.size() << "\t" << env::cell_size_x.size() << '\t' <<  env::neighbour_list_end_index.size() << std::endl;
       // This is massively inefficient for large numbers of cells

    //   std::cin.get();


       //calcualtes me

         //double m_e = pow((env::shield_Tc[0]-sim::temperature)/(env::shield_Tc[0]),0.365);


       //adds all cells which are not within the atomistic section to the the none atomsitic cells list or the atomistic cells list
       for (int cell = 0; cell < env::num_cells; cell++){


          // calculate if cell is part of shield structure
          bool included = true;//internal::in_shield(env::cell_coords_array_x[cell], env::cell_coords_array_y[cell], env::cell_coords_array_z[cell]);

             if (env::env_cell_is_in_atomistic_region[cell] == 0 && included){
                env::none_atomistic_cells.push_back(cell);
                env::num_env_cells ++;
          //      std::cout << cell << '\t' << env::num_env_cells << "\t" << env::none_atomistic_cells[env::num_env_cells -2]<< std::endl;
             }
             // if it is atomistic then add it to list of atomistic cells
             else if(env::env_cell_is_in_atomistic_region[cell] == 1){
                env::atomistic_cells.push_back(cell);
             }
             // otherwise don't add it to anything
       }
//
//       std::cout << "worked out which cells are atomistic/env" << "\t" << env::Ms.size() <<env::initial_spin_x.size() << "\t" <<std::endl;
       // This is massively inefficient for large numbers of cells

      // std::cin.get();

       //initialise the direction of the cell magnetisation
       //if the initial spin configuration is set to random - give each cell a random magnetisation

          for (int cl = 0; cl < env::num_env_cells; cl++){
             int cell = env::none_atomistic_cells[cl];
             int shield = env::shield_number[cell];
            // std::cout << shield << "\t" << cell << "\t" << env::random_spins.size() << "\t" << env::x_mag_array.size() << std::endl;
             // if (env::random_spins[shield]){
             // //if within the system
             //    env::x_mag_array[cell] = m_e*mtrandom::gaussian()*env::Ms[cell];
             //    env::y_mag_array[cell] = m_e*mtrandom::gaussian()*env::Ms[cell];
             //    env::z_mag_array[cell] = m_e*mtrandom::gaussian()*env::Ms[cell];
             //  }
             //  else{

                env::x_mag_array[cell] = env::initial_spin_x[shield]*env::Ms[cell];
                env::y_mag_array[cell] = env::initial_spin_y[shield]*env::Ms[cell];
                env::z_mag_array[cell] = env::initial_spin_z[shield]*env::Ms[cell];
            //    std::cout << env::x_mag_array[cell] << "\t" << env::Ms[cell] <<  std::endl;
             // }
          }
//
//
  //        std::cout << "set M" << "\t" << std::endl;//env::num_cells <<env::cell_coords_array_x.size() << "\t" << env::cell_size_x.size() << '\t' <<  env::neighbour_list_end_index.size() << std::endl;

       env::one_o_chi_para.resize(env::num_env_cells,0.1);
       env::one_o_chi_perp.resize(env::num_env_cells,0.1);

    //   std::cout << env::num_env_cells << "\t" << env::atomistic_cells.size() << '\t' << '\t' << env::num_cells << std::endl;

      //env::bias_shields();

    //initalise the demag fields
    env::initialise_demag_fields();

       // std::ofstream mfile;
       // mfile.open("m3.txt");
       //
       // for (int cell = 0; cell < cells::num_cells; cell++)
       // mfile<< cells::pos_and_mom_array[4*cell+0] + env::shift[0]<< '\t' << cells::pos_and_mom_array[4*cell+1]  +env::shift[1]<< '\t' << cells::pos_and_mom_array[4*cell+2]+env::shift[2] << '\t' << cells::mag_array_x[cell] <<'\t' << cells::mag_array_y[cell] <<'\t' << cells::mag_array_z[cell] <<std::endl;
       //
       //
       //
       //
       // std::ofstream pfile;
       // pfile.open("m2.txt");
       //
       // for (int cell = 0; cell < env::num_cells; cell++){
       //    pfile << cell << "\t" << env::cell_coords_array_x[cell]<< "\t" << env::cell_coords_array_y[cell]<< "\t" << env::cell_coords_array_z[cell]  <<  "\t" <<  env::x_mag_array[cell] << '\t' << env::y_mag_array[cell] << '\t' << env::z_mag_array[cell] <<std::endl;
       //  }
//          std::cin.get();


      return;

   }

} // end of environment namespace


void neighbour(int env, int mm, double overlap){
   environment::list_of_mm_cells_with_neighbours.push_back(mm);
   environment::list_of_env_cells_with_neighbours.push_back(env);
   environment::list_of_overlap_area.push_back(overlap);
   environment::num_interactions++;
   }
