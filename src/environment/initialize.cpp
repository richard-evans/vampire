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
#include "cells.hpp"
#include "vmpi.hpp"

#include <iostream>
#include <math.h>



std::vector <double> tmp_x;
std::vector <double> tmp_y;
std::vector <double> tmp_z;
std::vector <double> size_x;
std::vector <double> size_y;
std::vector <double> size_z;

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

       double mm_size_x, mm_size_y,mm_size_z;
       double env_size_x, env_size_y,env_size_z;

       if (env::shift[0] >= 0) {
          mm_size_x = system_dimensions_x + env::shift[0] +0.01;
          env_size_x = env::dim[0] + 0.01;
       }
       else{
          mm_size_x = system_dimensions_x +0.01;
          env_size_x = env::dim[0] + 0.01 + env::shift[0];
       }

       if (env::shift[1] >= 0) {
          mm_size_y = system_dimensions_y + env::shift[1] +0.01;
          env_size_y = env::dim[1] + 0.01;
       }
       else{
          mm_size_y = system_dimensions_y +0.01;
          env_size_y = env::dim[1] + 0.01 + env::shift[1];
       }


       if (env::shift[2] >= 0) {
          mm_size_z = system_dimensions_z + env::shift[2] +0.01;
          env_size_z = env::dim[2] + 0.01;
       }
       else{
          mm_size_z = system_dimensions_z +0.01;
          env_size_z = env::dim[2] + 0.01 + env::shift[2] ;
       }


       double dx,dy,dz;

       if (env_size_x < mm_size_x) dx = mm_size_x;
       else dx = env_size_x;
       if (env_size_y < mm_size_y) dy = mm_size_y;
       else dy = env_size_y;
       if (env_size_z < mm_size_z) dz = mm_size_z;
       else dz = env_size_z;
       //convert to the total system size

       // determine number of stacks in x and y (global)
       env::num_cells_x =  dx/env::cell_size;
       env::num_cells_y =  dy/env::cell_size;
       env::num_cells_z =  dz/env::cell_size;
       std::cout << env::num_cells_x << '\t' << env::num_cells_y << '\t' << env::num_cells_z <<std::endl;

       //total number of cells and cell volume
      // env::num_cells = env::num_cells_x*env::num_cells_y*env::num_cells_z;
       //env::cell_volume = env::cell_size*env::cell_size*env::cell_size;


       //stores the max and min coordinates of all env cells
       std::vector <double > x_max;
       std::vector <double > y_max;
       std::vector <double > z_max;
       std::vector <double > x_min;
       std::vector <double > y_min;
       std::vector <double > z_min;

       std::ofstream pfile2;
       pfile2.open("m2.txt");

      std::vector < std::vector < double> >corners;
      corners.resize(8);
      for (int i = 0; i < 8; ++i) corners[i].resize(3,0.0);

       //calculates the minimum and maximum positions of each cell for calcualtions of which cell the atomistic atoms are in later
       int cell = 0;

       for (int x = 0; x < env::num_cells_x; x++){
          for (int y = 0; y < env::num_cells_y; y++){
             for (int z = 0; z < env::num_cells_z; z++){

                double current_cell_size = env::cell_size;

                double pos_x = current_cell_size*x + current_cell_size/2.0;
                double pos_y = current_cell_size*y + current_cell_size/2.0;
                double pos_z = current_cell_size*z + current_cell_size/2.0;

                size_x.push_back(current_cell_size);
                size_y.push_back(current_cell_size);
                size_z.push_back(current_cell_size);
                tmp_x.push_back(pos_x);
                tmp_y.push_back(pos_y);
                tmp_z.push_back(pos_z);

              }
            }
          }
          int n_cell = 0;
          for (int cell = 0; cell < tmp_x.size(); cell ++ ){


                corners = env::calculate_corners(tmp_x[cell],tmp_y[cell],tmp_z[cell],size_x[cell],size_y[cell],size_z[cell]);
                int N = 0;
                int shield = 0;
                bool entered = false;
                for (int corner = 0; corner < 8; corner ++){

                  int within = env::in_shield(corners[corner][0] , corners[corner][1] , corners[corner][2]);
                  if (within != 0 && entered  == false){
                      N++;
                      shield = within;
                      entered = true;
                  }
                  else if (within != 0 && entered  == true){
                    if (within == shield){
                      N++;
                    }

                  }
                }

                if (N != 8 && size_x[cell] > 40.0){
            //      std::cout << "splitting cell \t:" << cell << '\t' <<  tmp_x[cell] << '\t' << tmp_y[cell] << '\t' << tmp_z[cell]  << '\t' << size_x[cell] << std::endl;
                  split_cells(tmp_x[cell],tmp_y[cell],tmp_z[cell], size_x[cell]);
            //      std::cout << "cell split" <<std::endl;
                }
                else if (N != 0 ) {

                  env::cell_coords_array_x.push_back(tmp_x[cell]);
                  env::cell_coords_array_y.push_back(tmp_y[cell]);
                  env::cell_coords_array_z.push_back(tmp_z[cell]);
                  env::cell_size_x.push_back(size_x[cell]);
                  env::cell_size_y.push_back(size_y[cell]);
                  env::cell_size_z.push_back(size_z[cell]);
                  env::cell_volume.push_back(size_x[cell]*size_y[cell]*size_z[cell]);
                  env::Ms.push_back(env::Ms_value/env::cell_volume[n_cell]);
                  env::ku.push_back(-env::ku_value/env::cell_volume[n_cell]);
                  x_max.push_back(tmp_x[cell] + size_x[cell]/2.0);
                  x_min.push_back(tmp_x[cell] - size_x[cell]/2.0);
                  y_max.push_back(tmp_y[cell] + size_y[cell]/2.0);
                  y_min.push_back(tmp_y[cell] - size_y[cell]/2.0);
                  z_max.push_back(tmp_z[cell] + size_z[cell]/2.0);
                  z_min.push_back(tmp_z[cell] - size_z[cell]/2.0);
                  pfile2 << N << '\t' << cell << '\t' <<  tmp_x[cell] << '\t' << tmp_y[cell] << '\t' << tmp_z[cell] << "\t" << size_x[cell] << '\t' << size_y[cell]<< '\t' << size_z[cell] <<  std::endl;
                  n_cell ++;

                }
              }


       env::num_cells = n_cell;
       std::cout << "Number of environment cells: " << env::num_cells <<  "\t cell size: " <<   env::cell_volume[0] << "\t" << env::cell_size << std::endl;

       //convert Ms from input to Ms = ms/V and Ku = ku/V



       //resize arrays
       env::x_mag_array.resize(env::num_cells,0.0);
       env::y_mag_array.resize(env::num_cells,0.0);
       env::z_mag_array.resize(env::num_cells,0.0);

       env::bias_field_x.resize(env::num_cells,0.0);
       env::bias_field_y.resize(env::num_cells,0.0);
       env::bias_field_z.resize(env::num_cells,0.0);

       env::dipole_field_x.resize(env::num_cells,0.0);
       env::dipole_field_y.resize(env::num_cells,0.0);
       env::dipole_field_z.resize(env::num_cells,0.0);

       env::neighbour_list_start_index.resize(env::num_cells,0.0);
       env::neighbour_list_end_index.resize(env::num_cells,0.0);

       env::list_env_cell_atomistic_cell.resize(cells::num_cells,0);
       env::env_cell_is_in_atomistic_region.resize(env::num_cells,0);

       environment_field_x.resize(cells::num_cells,0.0);
       environment_field_y.resize(cells::num_cells,0.0);
       environment_field_z.resize(cells::num_cells,0.0);



       std::cout << "Identifying environment cells which are atomistic" << std::endl;

       //loops over all atomistic cells to determine if the atomsitic simulation lies within an environment cell

       // std::ofstream pfile;
       // pfile.open("m3.txt");
       //
       // for(int lc=0; lc<cells::num_local_cells; lc++){
       //    int cell = cells::cell_id_array[lc];
       //    pfile << cell << "\t" << cells::pos_and_mom_array[4*cell + 0] << "\t" << cells::pos_and_mom_array[4*cell + 1] << "\t" << cells::pos_and_mom_array[4*cell + 2] <<  "\t" <<  cells::mag_array_x[cell] << '\t' << cells::mag_array_y[cell] << '\t' << cells::mag_array_z[cell] <<std::endl;
       //
       // }

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
            env::Ms.push_back(env::Ms_value/(cells::macro_cell_size_x*cells::macro_cell_size_y*cells::macro_cell_size_z));
            env::ku.push_back(-env::ku_value/(cells::macro_cell_size_x*cells::macro_cell_size_y*cells::macro_cell_size_z));
            env::list_env_cell_atomistic_cell[cell] = env_cell;
            //this cell is an atomistic cell.
            env::env_cell_is_in_atomistic_region.push_back(1);
            env::x_mag_array.push_back(cells::mag_array_x[cell]);
            env::y_mag_array.push_back(cells::mag_array_y[cell]);
            env::z_mag_array.push_back(cells::mag_array_z[cell]);
            pfile2 << 234 << '\t' << env_cell << '\t' <<  x << '\t' << y << '\t' << z << "\t" << cells::macro_cell_size_x << '\t' << cells::macro_cell_size_y << '\t' << cells::macro_cell_size_z <<  std::endl;
          }
       }

    //   std::cin.get();

        // std::ofstream pfile2;
        // pfile2.open("m2.txt");
       //
       // for (int cell = 0; cell < env::num_cells; cell++){
       //    pfile2 << cell << "\t" << env::cell_coords_array_x[cell] << "\t" << env::cell_coords_array_y[cell] << "\t" << env::cell_coords_array_z[cell] <<  "\t" <<  env::x_mag_array[cell] << '\t' << env::y_mag_array[cell] << '\t' << env::z_mag_array[cell] <<std::endl;
       //  }

       //calcualtes the neighbour lists for each cell.
       //if cells are neighbours add them to the neighbour list array_index
       //each cell has a start and end index index. the neighbours for each cell are within these index

       std::cout << "Calculating neighbour list for environment cells" << std::endl;
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
      //          std::cout << celli << '\t' << cellj << std::endl;
             }
          }
       }
    //   std::cin.get();


       //calcualtes me
       double m_e = pow((env::Tc-sim::temperature)/(env::Tc),0.365);

       //adds all cells which are not within the atomistic section to the the none atomsitic cells list or the atomistic cells list
       for (int cell = 0; cell < env::num_cells; cell++){

          // calculate if cell is part of shield structure
          bool included = internal::in_shield(env::cell_coords_array_x[cell], env::cell_coords_array_y[cell], env::cell_coords_array_z[cell]);
        //  std::cout << env::cell_coords_array_x[cell] << '\t' <<  env::cell_coords_array_y[cell] << '\t' <<  env::cell_coords_array_z[cell] << "\t" << included << std::endl;
          //std::cout << env::shift[0] << '\t' << env::cell_coords_array_y[cell] << '\t' << env::shift[2] << '\t' << env_size_x << '\t' << env_size_y << '\t' << env_size_z << '\t' << std::endl;
          // if (env::shift[0] < 0 && env::cell_coords_array_x[cell] < -env::shift[0]) included = false;
          // else if (env::shift[0] >= 0 && env::cell_coords_array_x[cell] > env_size_x + 15)  included = false;
          //
          // if (env::shift[1] >= 0 && env::cell_coords_array_y[cell] > env_size_y + 15)  included = false;
          // else if (env::shift[1] < 0 && env::cell_coords_array_y[cell] < -env::shift[1]) included = false;
          //
          // if (env::shift[2] >= 0 && env::cell_coords_array_z[cell] > env_size_z + 15)  included = false;
          // else if (env::shift[2] < 0 && env::cell_coords_array_z[cell] < -env::shift[2]) included = false;
             // check environment cell is in shields and not atomistic


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
      // std::cin.get();

       //initialise the direction of the cell magnetisation
       //if the initial spin configuration is set to random - give each cell a random magnetisation
       if (env::random_spins){

          for (int cl = 0; cl < env::num_env_cells; cl++){
             int cell = env::none_atomistic_cells[cl];
             //if within the system
                env::x_mag_array[cell] = m_e*mtrandom::gaussian()*env::Ms[cell];
                env::y_mag_array[cell] = m_e*mtrandom::gaussian()*env::Ms[cell];
                env::z_mag_array[cell] = m_e*mtrandom::gaussian()*env::Ms[cell];
          }
       }
       //if initial cell mangetisation is set as a direction set to that direction (normalised)
       else{
          //const double normal = sqrt(env::initial_spin[0]*env::initial_spin[0] + env::initial_spin[1]*env::initial_spin[1]  + env::initial_spin[2]*env::initial_spin[2]);

          for (int cl = 0; cl < env::num_env_cells; cl++){
             int cell = env::none_atomistic_cells[cl];
             //if within the system
             //std::cout << cell << env::x_mag_array[cell] << '\t' << env::y_mag_array[cell] << '\t' << env::z_mag_array[cell] <<std::endl;
               env::x_mag_array[cell] = env::initial_spin[0]*env::Ms[cell];
               env::y_mag_array[cell] = env::initial_spin[1]*env::Ms[cell];
               env::z_mag_array[cell] = env::initial_spin[2]*env::Ms[cell];

          }
       }

      // std::cout << env::num_env_cells << "\t" << env::atomistic_cells.size() << '\t' << '\t' << env::num_cells << std::endl;

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

void split_cells(double pos_x, double pos_y, double pos_z, double current_cell_size){

  for (int n = 0; n <8; n++ ){
    double new_cell_size = current_cell_size/2.0;
    size_x.push_back(new_cell_size);
    size_y.push_back(new_cell_size);
    size_z.push_back(new_cell_size);
    tmp_x.push_back(pos_x + x_cubes[n]*new_cell_size/2.0);
    tmp_y.push_back(pos_y + y_cubes[n]*new_cell_size/2.0);
    tmp_z.push_back(pos_z + z_cubes[n]*new_cell_size/2.0);
  //  std::cout << new_cell_size <<std::endl;
  }

}

void neighbour(int env, int mm, double overlap){
   environment::list_of_mm_cells_with_neighbours.push_back(mm);
   environment::list_of_env_cells_with_neighbours.push_back(env);
   environment::list_of_overlap_area.push_back(overlap);
   environment::num_interactions++;
   }
