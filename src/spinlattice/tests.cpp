//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Mara Strungaru 2022. All rights reserved.
//
//   Email: mara.strungaru@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
// Vampire headers
#include "sld.hpp"
#include "atoms.hpp"
#include "iostream"
#include "neighbours.hpp"
#include "material.hpp"
#include "constants.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>


// sld module headers
#include "internal.hpp"


namespace sld{

   //----------------------------------------------------------------------------
   // Function to initialize sld module
   //----------------------------------------------------------------------------
   void tests(){
      std::cout<<"*******************************************************"<<std::endl;
      std::cout<<"Testing the module of Spin-lattice dynamics simulations:"<<std::endl;
      std::cout<<"*******************************************************"<<std::endl;
      std::cout<<std::endl;
      std::cout<<std::endl;
      std::cout<<"*******************************************************"<<std::endl;
      std::cout<<"Testing forces and fields"<<std::endl;
      std::cout<<"*******************************************************"<<std::endl;

      //defining test quantites
      std::vector<double> potential_engt;
      potential_engt.resize(atoms::num_atoms,0);


      std::vector<double> xt_coord_array;
      std::vector<double> yt_coord_array;
      std::vector<double> zt_coord_array;

      xt_coord_array.resize(atoms::num_atoms,0.0);
      yt_coord_array.resize(atoms::num_atoms,0.0);
      zt_coord_array.resize(atoms::num_atoms,0.0);

      std::vector<double> xt_spin_array;
      std::vector<double> yt_spin_array;
      std::vector<double> zt_spin_array;

      xt_spin_array.resize(atoms::num_atoms,0.0);
      yt_spin_array.resize(atoms::num_atoms,0.0);
      zt_spin_array.resize(atoms::num_atoms,1.0);

      std::vector<double> forces_array_xt;
      std::vector<double> forces_array_yt;
      std::vector<double> forces_array_zt;

      forces_array_xt.resize(atoms::num_atoms,0.0);
      forces_array_yt.resize(atoms::num_atoms,0.0);
      forces_array_zt.resize(atoms::num_atoms,0.0);

      std::vector<double> fields_array_xt;
      std::vector<double> fields_array_yt;
      std::vector<double> fields_array_zt;

      fields_array_xt.resize(atoms::num_atoms,0.0);
      fields_array_yt.resize(atoms::num_atoms,0.0);
      fields_array_zt.resize(atoms::num_atoms,0.0);

      std::ofstream ofile_c("test_coupling.txt");
      std::ofstream ofile_e("test_exchange.txt");
      std::ofstream ofile_p("test_potential.txt");
      std::ofstream ofile_s("test_spins.txt");


      //set full printing precision
      std::cout<<std::setprecision(17)<<std::endl;

      int test_id;

      for( int i = 0; i < atoms::num_atoms; i++){

         xt_coord_array[i]=atoms::x_coord_array[i];
         yt_coord_array[i]=atoms::y_coord_array[i];
         zt_coord_array[i]=atoms::z_coord_array[i];

         if((xt_coord_array[i]-0.5*cs::system_dimensions[0]<1e-2)&&(yt_coord_array[i]-0.5*cs::system_dimensions[1]<1e-2)&&(zt_coord_array[i]-0.5*cs::system_dimensions[2]<1e-2)) test_id=i;

      }

      //choosing an atom and displace it inside the system
      std::cout<<"testing atom "<<test_id<<"\t"<<xt_coord_array[test_id]<<"\t"<<yt_coord_array[test_id]<<"\t"<<zt_coord_array[test_id]<<std::endl;
      for( int i = -1; i <= 1; i++){
         for( int j=-1; j <= 1; j++){
            //for (int k=-10; k<=10; k++){
               xt_coord_array[test_id]=atoms::x_coord_array[test_id]+i*0.1;
               yt_coord_array[test_id]=atoms::y_coord_array[test_id]+j*0.1;
               //zt_coord_array[test_id]=atoms::z_coord_array[test_id]+k*0.1;
               std::fill(fields_array_xt.begin(), fields_array_xt.begin()+atoms::num_atoms, 0.0);
               std::fill(fields_array_yt.begin(), fields_array_yt.begin()+atoms::num_atoms, 0.0);
               std::fill(fields_array_zt.begin(), fields_array_zt.begin()+atoms::num_atoms, 0.0);

               std::fill(forces_array_xt.begin(), forces_array_xt.begin()+atoms::num_atoms, 0.0);
               std::fill(forces_array_yt.begin(), forces_array_yt.begin()+atoms::num_atoms, 0.0);
               std::fill(forces_array_zt.begin(), forces_array_zt.begin()+atoms::num_atoms, 0.0);


               sld::internal::compute_forces_harmonic(test_id, test_id+1,
                                              atoms::neighbour_list_start_index, atoms::neighbour_list_end_index,
                                              atoms::type_array, atoms::neighbour_list_array,
                                              sld::internal::x0_coord_array, sld::internal::y0_coord_array, sld::internal::z0_coord_array,
                                              xt_coord_array, yt_coord_array, zt_coord_array,
                                              forces_array_xt, forces_array_yt, forces_array_zt, potential_engt);
               ofile_p<<std::setprecision(17)<<xt_coord_array[test_id]<<"\t"<<yt_coord_array[test_id]<<"\t"<<zt_coord_array[test_id]<<"\t"<<forces_array_xt[test_id]<<"\t"<<forces_array_yt[test_id]<<"\t"<<forces_array_zt[test_id]<<"\t"<<fields_array_xt[test_id]<<"\t"<<fields_array_yt[test_id]<<"\t"<<fields_array_zt[test_id]<<std::endl;

               std::fill(fields_array_xt.begin(), fields_array_xt.begin()+atoms::num_atoms, 0.0);
               std::fill(fields_array_yt.begin(), fields_array_yt.begin()+atoms::num_atoms, 0.0);
               std::fill(fields_array_zt.begin(), fields_array_zt.begin()+atoms::num_atoms, 0.0);

               std::fill(forces_array_xt.begin(), forces_array_xt.begin()+atoms::num_atoms, 0.0);
               std::fill(forces_array_yt.begin(), forces_array_yt.begin()+atoms::num_atoms, 0.0);
               std::fill(forces_array_zt.begin(), forces_array_zt.begin()+atoms::num_atoms, 0.0);


               sld::internal::compute_sld_coupling(test_id, test_id+1,
                                              atoms::neighbour_list_start_index, atoms::neighbour_list_end_index,
                                              atoms::type_array, atoms::neighbour_list_array,
                                              xt_coord_array, yt_coord_array, zt_coord_array,
                                              xt_spin_array, yt_spin_array, zt_spin_array,
                                              forces_array_xt, forces_array_yt, forces_array_zt,
                                              fields_array_xt, fields_array_yt, fields_array_zt);

               ofile_c<<std::setprecision(17)<<xt_coord_array[test_id]<<"\t"<<yt_coord_array[test_id]<<"\t"<<zt_coord_array[test_id]<<"\t"<<forces_array_xt[test_id]<<"\t"<<forces_array_yt[test_id]<<"\t"<<forces_array_zt[test_id]<<"\t"<<fields_array_xt[test_id]<<"\t"<<fields_array_yt[test_id]<<"\t"<<fields_array_zt[test_id]<<std::endl;

               std::fill(fields_array_xt.begin(), fields_array_xt.begin()+atoms::num_atoms, 0.0);
               std::fill(fields_array_yt.begin(), fields_array_yt.begin()+atoms::num_atoms, 0.0);
               std::fill(fields_array_zt.begin(), fields_array_zt.begin()+atoms::num_atoms, 0.0);

               std::fill(forces_array_xt.begin(), forces_array_xt.begin()+atoms::num_atoms, 0.0);
               std::fill(forces_array_yt.begin(), forces_array_yt.begin()+atoms::num_atoms, 0.0);
               std::fill(forces_array_zt.begin(), forces_array_zt.begin()+atoms::num_atoms, 0.0);


               sld::internal::compute_exchange(test_id, test_id+1,
                                              atoms::neighbour_list_start_index, atoms::neighbour_list_end_index,
                                              atoms::type_array, atoms::neighbour_list_array,
                                              xt_coord_array, yt_coord_array, zt_coord_array,
                                              xt_spin_array, yt_spin_array, zt_spin_array,
                                              forces_array_xt, forces_array_yt, forces_array_zt,
                                              fields_array_xt, fields_array_yt, fields_array_zt);

               ofile_e<<std::setprecision(17)<<xt_coord_array[test_id]<<"\t"<<yt_coord_array[test_id]<<"\t"<<zt_coord_array[test_id]<<"\t"<<forces_array_xt[test_id]<<"\t"<<forces_array_yt[test_id]<<"\t"<<forces_array_zt[test_id]<<"\t"<<fields_array_xt[test_id]<<"\t"<<fields_array_yt[test_id]<<"\t"<<fields_array_zt[test_id]<<std::endl;
               ofile_s<<std::setprecision(17)<<xt_coord_array[test_id]<<"\t"<<yt_coord_array[test_id]<<"\t"<<zt_coord_array[test_id]<<"\t"<<xt_spin_array[test_id]<<"\t"<<yt_spin_array[test_id]<<"\t"<<zt_spin_array[test_id]<<std::endl;

            }
         }
      //}


      ofile_c.close();
      ofile_p.close();
      ofile_e.close();
      ofile_s.close();

      return;


   }


} // end of sld namespace
