//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2016. All rights reserved.
//
//   Email: sj681@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>

// Vampire headers
#include "micromagnetic.hpp"
#include "atoms.hpp"
#include "material.hpp"
#include "cells.hpp"
#include "sim.hpp"
#include "../create/internal.hpp" // please fix
#include "vmpi.hpp"
#include "vio.hpp"
#include "create.hpp"

// micromagnetic module headers
#include "internal.hpp"

namespace mm = micromagnetic::internal;

namespace micromagnetic{

//----------------------------------------------------------------------------
// Function to initialize micromagnetic module
//----------------------------------------------------------------------------
void initialize(int num_local_cells,
                int num_cells,
                int num_atoms,
                int num_materials,
                std::vector<int> cell_array,                     //1D array storing which cell each atom is in
                std::vector<int> neighbour_list_array,           //1D vector listing the nearest neighbours of each atom
                std::vector<int> neighbour_list_start_index,     //1D vector storing the start index for each atom in the neighbour_list_array
                std::vector<int> neighbour_list_end_index,       //1D vector storing the end index for each atom in the neighbour_list_array
                const std::vector<int> type_array,               //1D array storing which material each cell is
                std::vector <mp::materials_t> material,          //1D vector of type material_t stiring the material properties
                std::vector <double> x_coord_array,
                std::vector <double> y_coord_array,
                std::vector <double> z_coord_array,
                std::vector <double> volume_array,               //1D vector storing the volume of each cell
                double Temperature,
                double num_atoms_in_unit_cell,
                double system_dimensions_x,
                double system_dimensions_y,
                double system_dimensions_z,
                std::vector<int> local_cell_array){


   if (micromagnetic::discretisation_type != 0){
   // Output informative message to user
   std::cout << "Initialising micromagnetic module" << std::endl;

   // initialise data output file
   mm::mm_output.open("mm_output.txt");

   // Determine number of local atoms for parallel computation
   int num_atoms_interactions;

   #ifdef MPICF
      num_atoms_interactions = vmpi::num_core_atoms+vmpi::num_bdry_atoms + vmpi::num_halo_atoms;
      num_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms; // only consider local atoms in MPI
   #else
      num_atoms_interactions = num_atoms;
   #endif

   mm::A.resize(num_cells*num_cells,0.0);
   mm::alpha.resize(num_cells,0.0);
   mm::one_o_chi_perp.resize(num_cells,0.0);
   mm::one_o_chi_para.resize(num_cells,0.0);
   mm::gamma.resize(num_cells,0.0);
   mm::ku.resize(num_cells,0.0);
   mm::ku_x.resize(num_cells,0.0);
   mm::ku_y.resize(num_cells,0.0);
   mm::ku_z.resize(num_cells,0.0);
   mm::ms.resize(num_cells,0.0);
   mm::T.resize(num_cells,0.0);
   mm::Tc.resize(num_cells,0.0);
   mm::alpha_para.resize(num_cells,0.0);
   mm::alpha_perp.resize(num_cells,0.0);
   mm::m_e.resize(num_cells,0.0);
   mm::macro_neighbour_list_start_index.resize(num_cells,0.0);
   mm::macro_neighbour_list_end_index.resize(num_cells,0.0);
   micromagnetic::cell_discretisation_micromagnetic.resize(num_cells,true);
   mm::ext_field.resize(3,0.0);
   mm::pinning_field_x.resize(num_cells,0.0);
   mm::pinning_field_y.resize(num_cells,0.0);
   mm::pinning_field_z.resize(num_cells,0.0);
   mm::cell_material_array.resize(num_cells,0.0);


   // These functions vectors with the parameters calcualted from the function
   mm::ms =                   mm::calculate_ms(num_local_cells,num_atoms,num_cells, cell_array, type_array,material,local_cell_array);
   mm::alpha =                mm::calculate_alpha(num_local_cells,num_atoms, num_cells, cell_array, type_array, material,local_cell_array);
   mm::Tc =                   mm::calculate_tc(num_local_cells, local_cell_array,num_atoms_interactions, num_cells, cell_array,neighbour_list_array,
                                               neighbour_list_start_index, neighbour_list_end_index, type_array, material);
   mm::ku =                   mm::calculate_ku(num_atoms, num_cells, cell_array, type_array, material);
   mm::gamma =                mm::calculate_gamma(num_atoms, num_cells, cell_array,type_array,material,num_local_cells,local_cell_array);
   mm::A =                    mm::calculate_a(num_atoms_interactions, num_cells, num_local_cells,cell_array, neighbour_list_array, neighbour_list_start_index,
                                              neighbour_list_end_index, type_array,  material, volume_array, x_coord_array,
                                              y_coord_array, z_coord_array, num_atoms_in_unit_cell, local_cell_array);

   // calculate spin transfer torque parameters
   mm::calculate_stt(num_atoms, num_cells, cell_array, type_array, material, mm::stt_rj, mm::stt_pj);

std::ofstream ofile("initial_parameters.txt");

for (int lc = 0; lc < num_local_cells; lc++){
 int cell = local_cell_array[lc];
 ofile<<cell << '\t' <<  mm::ms[cell] << '\t' << mm::alpha[cell] << '\t' << mm::Tc[cell] << '\t' << mm::ku[cell] << '\t' << mm::gamma[cell] << std::endl;

}
for (int proc = 0; proc < vmpi::num_processors; proc++ ){
 if (vmpi::my_rank == proc)   std::cerr << proc << "\t" << mm::macro_neighbour_list_array.size() << "\t" <<  std::endl;

}




   //---------------------------------------------------------------------------
   // Loop over all cells to determine if a multiscale solver is required
   //---------------------------------------------------------------------------
   if (discretisation_type == 1){
      for (int lc = 0; lc < num_local_cells; lc++){
         int cell = local_cell_array[lc];
       //  std::cout <<type_array[cell] << "\t" << x_coord_array[cell] << '\t' <<y_coord_array[cell] << '\t' <<z_coord_array[cell] << '\t' <<  mm::ms[cell] << std::endl;
         if (mm::Tc[cell] < 0) {
            discretisation_type = 2;
         }
      }
   }

   //----------------------------------------------------------------------------------
   // if multiscale simulation work out which cells/atoms are micromagnetic/atomistic
   //
   // STILL NEED TO FIX THIS IN PARALLEL!!
   // Need to make sure that actual micromagnetic cells are correctly allocated to
   // processors
   //
   //----------------------------------------------------------------------------------
   if (discretisation_type == 2){
      //loops over all atoms and if any atom in the cell is atomistic the whole cell becomes atomistic else the cell is micromagnetic
      // this doesnt do what you think it does... last atom defines if mm/atomistic, not any...
      for (int atom =0; atom < num_atoms; atom++){
         int cell = cell_array[atom];
         int mat  = type_array[atom];
         micromagnetic::cell_discretisation_micromagnetic[cell] = mp::material[mat].micromagnetic_enabled;
         //unless the cell contains AFM atoms, then it is always atomistic
         if (mm::Tc[cell] < 0) micromagnetic::cell_discretisation_micromagnetic[cell] = 0;
      }

      //loops over all atoms saves each atom at micromagnetic or atomistic depending on whether the cell is microamgnetic or atomistic
      for (int atom =0; atom < num_atoms; atom++){
         int cell = cell_array[atom];
         //id atomistic add to numner of atomisic atoms
         if (micromagnetic::cell_discretisation_micromagnetic[cell] == 0) {
            list_of_atomistic_atoms.push_back(atom);
            number_of_atomistic_atoms++;
         }
         //if micromagnetic add to the numebr of micromagnetic cells.
         else {
            list_of_none_atomistic_atoms.push_back(atom);
            number_of_none_atomistic_atoms++;
         }
      }

      //if simulation is micromagnetic all cells are made micromagnetic cells
      for (int lc = 0; lc < num_local_cells; lc++){
         int cell = local_cell_array[lc];
         if (micromagnetic::cell_discretisation_micromagnetic[cell] == 1 && mm::ms[cell] > 1e-30) {
            list_of_micromagnetic_cells.push_back(cell);
            number_of_micromagnetic_cells ++;
         }
         // Otherwise list cells as non-magnetic (empty)
         else{
            list_of_empty_micromagnetic_cells.push_back(cell);
         }
      }

      // Need to allocate cells to processors

   }

   //-------------------------------------------------------------------------------------------
   // if micromagnetic simulation all cells are micromagnetic and all atoms are micromagnetic
   //-------------------------------------------------------------------------------------------
   else {

      // wait for all processors
      vmpi::barrier();

      // array to store list of magnetic cells
      std::vector<int> list_of_magnetic_cells(0);

      // determine total number of magnetic cells
      int num_magnetic_cells = 0;
      for(int c = 0; c < cells::num_cells; c++){
         if(mm::ms[c] > 1e-40 && mm::Tc[c] > 0.1){
            list_of_magnetic_cells.push_back(c);
            num_magnetic_cells++;
         }
         else{
            list_of_empty_micromagnetic_cells.push_back(c);
         }
      }

      // Output informative message to user
      std::cout << "Number of simulated micromagnetic cells: " << num_magnetic_cells << std::endl;
      zlog << zTs() << "Number of simulated micromagnetic cells: " << num_magnetic_cells << std::endl;

      // now divide cells over processors, allocating extra to the last processor
      int num_cells_pp = num_magnetic_cells / vmpi::num_processors;
      int my_num_cells = num_cells_pp;
      if(vmpi::my_rank == vmpi::num_processors - 1) my_num_cells = num_magnetic_cells - (vmpi::num_processors-1)*num_cells_pp;

      // loop over all magnetic cells, allocating cells to processors in linear fashion
      const int start = vmpi::my_rank*num_cells_pp;
      const int end   = vmpi::my_rank*num_cells_pp + my_num_cells;
      for(int c = start; c < end; c++){
         list_of_micromagnetic_cells.push_back( list_of_magnetic_cells[c] );
         number_of_micromagnetic_cells ++;
      }

      vmpi::barrier();

      for (int atom =0; atom < num_atoms; atom++){
         list_of_none_atomistic_atoms.push_back(atom);
         number_of_none_atomistic_atoms++;
      }

   }

// std::ofstream ofile2("exchange.txt");
// for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){
//   int cell = list_of_micromagnetic_cells[lc];
//    //loops over all other cells with interactions to this cell
//    const int start = mm::macro_neighbour_list_start_index[cell];
//    const int end = mm::macro_neighbour_list_end_index[cell] +1;
//    std::cerr << vmpi::my_rank << '\t' << cell << '\t' << start << '\t' << end << "\t" << number_of_micromagnetic_cells << std::endl;
// if (vmpi::my_rank ==0){
//    for(int j = start;j< end;j++){
//       const int cellj = mm::macro_neighbour_list_array[j];
//       // calculate reduced exchange constant factor
//      // ofile2 << cell << '\t' << cellj << '\t' << mm::A[j] <<std::endl;

//         }
//      }
// }

   //-------------------------------------------------------------------------------------------
   // for field calculations you need to access the atoms in numerically consecutive lists.
   // therefore you need to create lists of consecutive lists
   // loops over all atoms if atom is not one minus the previous atom then create a new list.
   //-------------------------------------------------------------------------------------------
   if (number_of_atomistic_atoms > 0){
      int end = list_of_atomistic_atoms[0];
      int begin = list_of_atomistic_atoms[0];
      for(int atom_list=1;atom_list<number_of_atomistic_atoms;atom_list++){
         int atom = list_of_atomistic_atoms[atom_list];
         int last_atom = list_of_atomistic_atoms[atom_list - 1];
         if ((atom != last_atom +1) || (atom_list == number_of_atomistic_atoms -1)){
            end = atom +1;
            mm::fields_neighbouring_atoms_begin.push_back(begin);
            mm::fields_neighbouring_atoms_end.push_back(end);
            begin = atom + 1;
         }
      }
   }

   //--------------------------------------------------------------------------------------------------
   // Pre-calculate micromagnetic parameters (not technically needed)
   //--------------------------------------------------------------------------------------------------
   //mm::calculate_chi_para(number_of_micromagnetic_cells, list_of_micromagnetic_cells, mm::one_o_chi_para, mm::T, mm::Tc);
   //mm::calculate_chi_perp(number_of_micromagnetic_cells, list_of_micromagnetic_cells, mm::one_o_chi_perp, mm::T, mm::Tc);

   //-------------------------------------------------------------------------------------------
   // Save the cell material type to enable setting material specific properties
   // taking average over constituent atoms
   //-------------------------------------------------------------------------------------------
   std::vector< std::vector < int> > counter;
   counter.resize(num_cells);
   for (int i = 0; i < num_cells; ++i)
    counter[i].resize(num_materials);

   for (int atom = 0; atom < num_atoms_interactions; atom ++){
      int mat = type_array[atom];
      int cell = cell_array[atom];
      counter[cell][mat] ++;
      // if ( mm::cell_material_array[cell] < mat){
      //   mm::cell_material_array[cell] = mat;
      // }
   }
   for (int cell = 0; cell < num_cells; cell ++ ){
            int largest = 0;
   for (int mat = 0; mat < num_materials; mat++){

      if (counter[cell][mat] > largest){
         mm::cell_material_array[cell] = mat;
         largest = counter[cell][mat];
      }
   }
   }
    // for (int cell = 0; cell < num_cells; cell++ ){
    //    int mat = mm::cell_material_array[cell];
    //    if (vmpi::my_rank == 0) std::cerr << cell << '\t' << mat << std::endl;
    //  }

   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &mm::cell_material_array[0],     num_cells,    MPI_DOUBLE,    MPI_MAX, MPI_COMM_WORLD);
   #endif

   // for (int cell = 0; cell < num_cells; cell++ ){
   //    int mat = mm::cell_material_array[cell];
   //    if (vmpi::my_rank == 1) std::cerr << cell << '\t' << mat << std::endl;
   //  }

   //--------------------------------------------------------------------------------------------------
   // Pre-calculate SAF properties
   //--------------------------------------------------------------------------------------------------

   mm::mat_vol.resize(mp::num_materials,0);
   mm::mat_ms.resize(mp::num_materials,0);
   mm::prefactor.resize(mp::num_materials,0);
   for (int cell = 0; cell < num_cells; cell++ ){
      int mat = mm::cell_material_array[cell];
      mm::mat_vol[mat] = mm::mat_vol[mat] + volume_array[cell];
      mm::mat_ms[mat] = mm::mat_ms[mat] + mm::ms[cell];
      // std::cout <<mm::cell_material_array[cell] << "\t" << x_coord_array[cell] << '\t' <<y_coord_array[cell] << '\t' <<z_coord_array[cell] << '\t' <<  mm::ms[cell] << std::endl;
      //std::cout <<
   }
   for (int mat = 0; mat < mp::num_materials; mat++ ){
      double min=create::internal::mp[mat].min;
      double max=create::internal::mp[mat].max;
      double t = max - min;
      t = t*system_dimensions_z;
      mm::prefactor[mat] = 1/(t*mm::mat_ms[mat]/mm::mat_vol[mat]);
   }

   //--------------------------------------------------------------------------------------------------
   //Replace atomsitic fields with SAF accross the boundary
   //--------------------------------------------------------------------------------------------------
   double area = cells::macro_cell_size_x*cells::macro_cell_size_y;
   for (int cell = 0; cell < num_cells; cell++ ){

      //double zi = cells::pos_and_mom_array[4*cell+2];
      const int mat = mm::cell_material_array[cell];
      const int start = mm::macro_neighbour_list_start_index[cell]; // save start index for neighbour list
      const int end = mm::macro_neighbour_list_end_index[cell] +1;  // save end index for neighbour list
         // loop over neighbouring cells
         for(int j = start;j< end;j++){

            const int cellj = mm::macro_neighbour_list_array[j];
            const int matj = mm::cell_material_array[cellj];
            // Check if spaced SAF is included


            if (mp::material[mat].enable_SAF && mp::material[matj].enable_SAF && mat != matj){
               // check that materials are different
            //std::cout <<"SAF" << mat << "\t" << matj << "\t" << mp::material[mat].enable_SAF << "\t" << mp::material[matj].enable_SAF << std::endl;
             //  std::cout << "enter2" << std::endl;
               double dx = cells::pos_array[cell*3 +0] - cells::pos_array[cellj*3 +0];
               double dy = cells::pos_array[cell*3 +1] - cells::pos_array[cellj*3 +1];
               //double dz = cells::pos_array[cell*3 +2] - cells::pos_array[cellj*3 +2];
               //if (mat == mm::resistance_layer_1 && matj == mm::resistance_layer_2 && dx*dx < cs::unit_cell.dimensions[0]*cs::unit_cell.dimensions[0] && dy*dy < cs::unit_cell.dimensions[1]*cs::unit_cell.dimensions[1]){
               if (dx*dx < cs::unit_cell.dimensions[0]*cs::unit_cell.dimensions[0] && dy*dy < cs::unit_cell.dimensions[1]*cs::unit_cell.dimensions[1]){

                  // why mj^1.66?? need to check how this actually works. why / ms[cell]?
                  //Ac = -prefactor[matj]*mp::material[mat].SAF[matj];
                  mm::A[j] = -area*mp::material[mat].SAF[matj]/mm::ms[cell];
                  //if (mm_correction == true) Ac = 2*Ac/cells::macro_cell_size[2];

               }
               else {
                  mm::A[j] = 0.0;
               }
            }


            // option to override atomistic exchange with a micromagnetic value
            if (mp::material[mat].override_atomsitic[matj] == true){
               //double Area = cells::macro_cell_size*cells::macro_cell_size;
               //double Volume = cells::macro_cell_size*cells::macro_cell_size*cells::macro_cell_size;
               //Ac = -2*pow(mj,1.66)*mp::material[mat].EF_MM[matj]/(ms[cell]*Area);
               mm::A[j] = 2.0*mp::material[mat].EF_MM[matj]/(mm::ms[cell]);
            }
            // Output SAF coupling for testing
            //if (mp::material[mat].enable_SAF && mp::material[matj].enable_SAF && mat != matj) std::cout << cell << '\t' << cellj << "\t" << mm::A[j]  <<std::endl;
         }
   }

   //--------------------------------------------------------------------------------------------------
   // Initialise pinning field calculation
   //--------------------------------------------------------------------------------------------------
   for (int cell = 0; cell < num_cells; cell++ ){

      //double zi = cells::pos_and_mom_array[4*cell+2];
      int mat = mm::cell_material_array[cell];

      if (mp::material[mat].pinning_field_unit_vector[0]+ mp::material[mat].pinning_field_unit_vector[1] + mp::material[mat].pinning_field_unit_vector[2]!= 0.0){
         //double Area = cells::macro_cell_size_x*cells::macro_cell_size_y;
         //  std::cout << mp::material[mat].pinning_field_unit_vector[0] << '\t' <<mp::material[mat].pinning_field_unit_vector[1] << '\t' << mp::material[mat].pinning_field_unit_vector[2] << '\t' << mm::ms[cell] << '\t' << Area << std::endl;
         mm::pinning_field_x[cell] = mm::prefactor[mat]*mp::material[mat].pinning_field_unit_vector[0];
         mm::pinning_field_y[cell] = mm::prefactor[mat]*mp::material[mat].pinning_field_unit_vector[1];
         mm::pinning_field_z[cell] = mm::prefactor[mat]*mp::material[mat].pinning_field_unit_vector[2];

     // std::cout << mm::prefactor[mat] << '\t' << mm::pinning_field_x[cell] << "\t" << mm::pinning_field_y[cell] << "\t" << mm::pinning_field_z[cell] << "\t" << std::endl;
         // std::cout <<prefactor*mp::material[mat].pinning_field_unit_vector[2] << '\t' << prefactor2*mp::material[mat].pinning_field_unit_vector[2] << '\t' << prefactor3*mp::material[mat].pinning_field_unit_vector[2] << '\t' << prefactor4*mp::material[mat].pinning_field_unit_vector[2] <<  std::endl;
         //               n_cells

      }
   }

   //--------------------------------------------------------------------------------------------------
   // Pinning field corrections
   //--------------------------------------------------------------------------------------------------
   if (mm::mm_correction == true){
      for (int cell = 0; cell < num_cells; cell++ ){
         mm::pinning_field_x[cell] = 2.0*mm::pinning_field_x[cell]/cells::macro_cell_size_x;
         mm::pinning_field_y[cell] = 2.0*mm::pinning_field_y[cell]/cells::macro_cell_size_y;
         mm::pinning_field_z[cell] = 2.0*mm::pinning_field_z[cell]/cells::macro_cell_size_z;
      }
   }
   // loop over all cells
  //  for (int cell = 0; cell < cells::num_cells; cell++ ){
  //          int mat = mm::cell_material_array[cell];
  // //   std::cout << cell << "\t" << mat << '\t' << mm::pinning_field_x[cell] << '\t' << mm::pinning_field_y[cell] << '\t' <<  mm::pinning_field_z[cell] << std::endl;
  //  }

   //--------------------------------------------------------------------------------------------------
   // Initialise restistance calculation
   //--------------------------------------------------------------------------------------------------
   //std::cout << "AAAJKALJK" << std::endl;
   if (enable_resistance && mm::resistance_layer_2 != mm::resistance_layer_1 ){

      // loop over all cells
     // std::cout <<"CELLS" <<  cells::num_cells << '\t' << cells::pos_and_mom_array.size() <<std::endl;
      for (int cell = 0; cell < cells::num_cells; cell++ ){

         int mat = mm::cell_material_array[cell];
         const int start = mm::macro_neighbour_list_start_index[cell];
         const int end = mm::macro_neighbour_list_end_index[cell] +1;
       //       std::cout << cells::pos_and_mom_array[cell*4 +0] << '\t' << cells::pos_and_mom_array[cell*4 +1] << '\t' << cells::pos_and_mom_array[cell*4 +2] <<  "\t" << cells::pos_and_mom_array[cell*4 +3] << "\t" << mat << std::endl;//'\t' <<cells::pos_and_mom_array[cellj*4 +0] << '\t' <<  cells::pos_and_mom_array[cellj*4 +1] << '\t' << cells::pos_and_mom_array[cellj*4 +2] << '\t' << std::endl;

         for(int j = start;j< end;j++){
            // calculate reduced exchange constant factor
            const int cellj = mm::macro_neighbour_list_array[j];
            int matj =mm::cell_material_array[cellj];
            double dx = cells::pos_array[cell*3 +0] - cells::pos_array[cellj*3 +0];
            double dy = cells::pos_array[cell*3 +1] - cells::pos_array[cellj*3 +1];
            //double dz = cells::pos_array[cell*3 +2] - cells::pos_array[cellj*3 +2];
//std::cout << cell << '\t' << cellj << '\t' << mat << '\t' << matj << "\t" << mm::resistance_layer_2 << '\t' << mm::resistance_layer_1 <<"\t" << cells::pos_array[cell*3 +0] << '\t' << cells::pos_array[cell*3 +1] << '\t' << cells::pos_array[cell*3 +2] << '\t' <<  cells::pos_array[cellj*3 +1] << '\t' << cells::pos_array[cellj*3 +2] << '\t' << dx << "\t" << dy << '\t' << dz << std::endl;
            if (mat == mm::resistance_layer_1 && matj == mm::resistance_layer_2 && dx*dx < cs::unit_cell.dimensions[0]*cs::unit_cell.dimensions[0] && dy*dy < cs::unit_cell.dimensions[1]*cs::unit_cell.dimensions[1]){
     //          std::cout << cells::pos_and_mom_array[cell*4 +0] << '\t' << cells::pos_and_mom_array[cell*4 +1] << '\t' << cells::pos_and_mom_array[cell*4 +2] << '\t' <<cells::pos_and_mom_array[cellj*4 +0] << '\t' <<  cells::pos_and_mom_array[cellj*4 +1] << '\t' << cells::pos_and_mom_array[cellj*4 +2] << '\t' << std::endl;
      //    std::cout <<  dx << '\t' << dy << "\t" << dz << "\t" << mat << "\t" << matj << std::endl;
            //         std::cout << mm::resistance_layer_1 << '\t' << mm::resistance_layer_2 <<std::endl;
        //       std::cout << cell << '\t' << cellj << "\t" <<mat << '\t' << matj << std::endl;// x_coord_array[cell] << "\t" <<y_coord_array[cell] << "\t" <<z_coord_array[cell] << "\t" <<  x_coord_array[cellj] << "\t" <<y_coord_array[cellj] << "\t" <<z_coord_array[cellj] << "\t" <<std::endl;
               mm::overlap_area = mm::overlap_area + cells::macro_cell_size_x*cells::macro_cell_size_y;
               std::cout << cell << '\t' << cellj << "\t" <<mm::overlap_area << std::endl;// x_coord_array[cell] << "\t" <<y_coord_array[cell] << "\t" <<z_coord_array[cell] << "\t" <<  x_coord_array[cellj] << "\t" <<y_coord_array[cellj] << "\t" <<z_coord_array[cellj] << "\t" <<std::endl;

            }
         }
      }
   }
   }
   //--------------------------------------------------------------------------------------------------
   // Initialise bias magnets
   //--------------------------------------------------------------------------------------------------
   if(mm::bias_magnets == true){
      mm::bias_field_x.resize(num_cells,0.0);
      mm::bias_field_y.resize(num_cells,0.0);
      mm::bias_field_z.resize(num_cells,0.0);
      atomistic_bias_field_x.resize(num_atoms,0.0);
      atomistic_bias_field_y.resize(num_atoms,0.0);
      atomistic_bias_field_z.resize(num_atoms,0.0);
      mm::calculate_bias_magnets(system_dimensions_x,system_dimensions_y,system_dimensions_z);
      //std::cin.get();
      // For MPI version, only add local atoms
       #ifdef MPICF
          int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
       #else
          int num_local_atoms = atoms::num_atoms;
       #endif

      for (int atom =0; atom < num_local_atoms; atom++){
         int cell = cell_array[atom];
         //std::cout << atom << '\t' << cell << '\t' <<  mm::bias_field_x[cell] << '\t' << atomistic_bias_field_x[atom] << std::endl;
         atomistic_bias_field_x[atom] = mm::bias_field_x[cell];
         atomistic_bias_field_y[atom] = mm::bias_field_y[cell];
         atomistic_bias_field_z[atom] = mm::bias_field_z[cell];
      //   std::cout << atom << '\t' << cell << '\t' <<  mm::bias_field_x[cell] << '\t' << atomistic_bias_field_x[atom] << std::endl;
      }

   }
   if (discretisation_type == 1){
      for(int i=0; i<cells::num_cells; ++i) {
         cells::mag_array_x[i] = 0.0;
         cells::mag_array_y[i] = 0.0;
         cells::mag_array_z[i] = 0.0;
      }

      #ifdef MPICF
         int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
      #else
         int num_local_atoms = num_atoms;
      #endif

      // calulate total moment in each cell
      for(int i=0;i<num_local_atoms;++i) {
         int cell = cells::atom_cell_id_array[i];
         int type = type_array[i];
         //// Consider only cells with n_atoms != 0
         //if(cells::num_atoms_in_cell[cell]>0){
            const double mus = mp::material[type].mu_s_SI;
            // Consider only magnetic elements
            if(mp::material[type].non_magnetic==0){
               cells::mag_array_x[cell] += atoms::x_spin_array[i]*mus;
               cells::mag_array_y[cell] += atoms::y_spin_array[i]*mus;
               cells::mag_array_z[cell] += atoms::z_spin_array[i]*mus;
            }
         //}
      }

      #ifdef MPICF
      // Reduce magnetisation on all nodes
      MPI_Allreduce(MPI_IN_PLACE, &cells::mag_array_x[0],   cells::mag_array_x.size(),   MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &cells::mag_array_y[0],   cells::mag_array_y.size(),   MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &cells::mag_array_z[0],   cells::mag_array_z.size(),   MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
      #endif
      }
   //}
    //  std::cin.get();
   return;

}

} // end of micromagnetic namespace
