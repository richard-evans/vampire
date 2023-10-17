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

//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Mara Stungaru 2016. All rights reserved.
//
//   Email: mss555@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "errors.hpp"
#include "exchange.hpp"
#include "unitcell.hpp"
#include "vio.hpp"
#include "atoms.hpp"
#include "create.hpp"
#include "material.hpp"


// exchange module headers
#include "internal.hpp"

// Vampire headers

namespace exchange{

namespace internal{

   void initialize_four_spin_exchange(std::vector<std::vector <neighbours::neighbour_t> >& cneighbourlist){


      // if four spin exchange is not needed then do nothing
      if(!internal::enable_fourspin) return;



      exchange::internal::four_spin_neighbour_list_start_index.resize(atoms::num_atoms,0);
      exchange::internal::four_spin_neighbour_list_end_index.resize(atoms::num_atoms,0);

      //distances here are for a cubic system of normalised dimension
      double nn_distance = internal::fs_cutoff_1;
      double nnn_distance = internal::fs_cutoff_2;
     
      double d1,d2,d3;
      double ucx,ucy,ucz; 
      double x_a,y_a,z_a;
      double x_b,y_b,z_b;
      double x_c,y_c,z_c;
      int k1=1, k2=1,k3=1;
      int counter=0;
      int counter_sort=0;
      int c1,c2;
      
      //vectors to store the nearest and next nearest neighbours
      std::vector <int> first_neigh(0);
      std::vector <int> start_first_neigh(0);
      std::vector <int> end_first_neigh(0);

      first_neigh.resize(20*atoms::num_atoms);
      start_first_neigh.resize(atoms::num_atoms);
      end_first_neigh.resize(atoms::num_atoms);

     //to print out the four-spin interaction
      std::ofstream ofile;
      ofile.open("fourspin_quartets.txt");
     
    
     ucx=1.0;
     ucy=1.0;
     ucz=1.0;




	int counter1=0;

    ///this part computes the first order exchange neighbourlist and second



    start_first_neigh[0]=0;


    for(int i=0;i< atoms::num_atoms ;i++){
       const int start=atoms::neighbour_list_start_index[i];
	   const int end=atoms::neighbour_list_end_index[i]+1;
	   int am=0;
	   double dist=0;	
	

     for(int a=start;a<end;a++){

         dist=std::sqrt(std::pow(((cneighbourlist[i][am].vx)/ucx),2)+std::pow(((cneighbourlist[i][am].vy)/ucy),2)+std::pow(((cneighbourlist[i][am].vz)/ucz),2));

         if ((dist <= nn_distance+0.01) && (dist >= nn_distance-0.01)){
             first_neigh[counter1]=cneighbourlist[i][am].nn;
             counter1=counter1+1;
            }


          am=am+1;


      }
      end_first_neigh[i]=counter1-1;
      start_first_neigh[i+1]=end_first_neigh[i]+1;

  }
	
	

    //now create the fourspin interactions
    int n_interactions = 0;


    for(int i=0;i< atoms::num_atoms ;i++){
        const int start=start_first_neigh[i];
	    const int end=end_first_neigh[i]+1;	
        const int imaterial = atoms::type_array[i];
        const double imus = 1.0 / mp::material[imaterial].mu_s_SI; // get inverse spin moment
        
	     four_spin_neighbour_list_start_index[i] = n_interactions;

	for(int a=start;a<end;a++){
	const int jmaterial = atoms::type_array[a];
	    for(int b=a;b<end;b++){
         	for(int c=b;c<end;c++){
         	
         	
                 //if system is in the first or second half of the system; this is for pbc (-1)**k1
                 if(atoms::x_coord_array[i] <= (cs::system_dimensions[0]/2.0)) k1=1; else k1=2;
                 if(atoms::y_coord_array[i] <= (cs::system_dimensions[1]/2.0)) k2=1; else k2=2;
                 if(atoms::z_coord_array[i] <= (cs::system_dimensions[2]/2.0)) k3=1; else k3=2;
                 
                 

  
                 //distances with pbc, this can be changed to be more simple
                 if (abs(atoms::x_coord_array[i]-atoms::x_coord_array[first_neigh[a]])/ucx >= 1.01) x_a=atoms::x_coord_array[first_neigh[a]]+std::pow(-1,k1)*cs::system_dimensions[0];
                 else x_a=atoms::x_coord_array[first_neigh[a]];

                 if (abs(atoms::y_coord_array[i]-atoms::y_coord_array[first_neigh[a]])/ucy >= 1.01) y_a=atoms::y_coord_array[first_neigh[a]]+std::pow(-1,k2)*cs::system_dimensions[1];
                 else y_a=atoms::y_coord_array[first_neigh[a]];
                
                 if (abs(atoms::z_coord_array[i]-atoms::z_coord_array[first_neigh[a]])/ucz >= 1.01) z_a=atoms::z_coord_array[first_neigh[a]]+std::pow(-1,k3)*cs::system_dimensions[2];
                 else z_a=atoms::z_coord_array[first_neigh[a]];

                 if (abs(atoms::x_coord_array[i]-atoms::x_coord_array[first_neigh[b]])/ucx >= 1.01) x_b=atoms::x_coord_array[first_neigh[b]]+std::pow(-1,k1)*cs::system_dimensions[0];
                 else x_b=atoms::x_coord_array[first_neigh[b]];
                 
                 if (abs(atoms::y_coord_array[i]-atoms::y_coord_array[first_neigh[b]])/ucy >= 1.01) y_b=atoms::y_coord_array[first_neigh[b]]+std::pow(-1,k2)*cs::system_dimensions[1];
                 else y_b=atoms::y_coord_array[first_neigh[b]];
                 
                 if (abs(atoms::z_coord_array[i]-atoms::z_coord_array[first_neigh[b]])/ucz >= 1.01) z_b=atoms::z_coord_array[first_neigh[b]]+std::pow(-1,k3)*cs::system_dimensions[2];
                 else z_b=atoms::z_coord_array[first_neigh[b]];

                 if (abs(atoms::x_coord_array[i]-atoms::x_coord_array[first_neigh[c]])/ucx >= 1.01) x_c=atoms::x_coord_array[first_neigh[c]]+std::pow(-1,k1)*cs::system_dimensions[0];
                 else x_c=atoms::x_coord_array[first_neigh[c]];
                 
                 if (abs(atoms::y_coord_array[i]-atoms::y_coord_array[first_neigh[c]])/ucy >= 1.01) y_c=atoms::y_coord_array[first_neigh[c]]+std::pow(-1,k2)*cs::system_dimensions[1];
                 else y_c=atoms::y_coord_array[first_neigh[c]];
                 
                 if (abs(atoms::z_coord_array[i]-atoms::z_coord_array[first_neigh[c]])/ucz >= 1.01) z_c=atoms::z_coord_array[first_neigh[c]]+std::pow(-1,k3)*cs::system_dimensions[2];
                 else z_c=atoms::z_coord_array[first_neigh[c]];

                 d1=std::sqrt(std::pow(((x_a-x_b)/ucx),2)+std::pow(((y_a-y_b)/ucy),2)+std::pow(((z_a-z_b)/ucz),2));
                 d2=std::sqrt(std::pow(((x_c-x_b)/ucx),2)+std::pow(((y_c-y_b)/ucy),2)+std::pow(((z_c-z_b)/ucz),2));
                 d3=std::sqrt(std::pow(((x_a-x_c)/ucx),2)+std::pow(((y_a-y_c)/ucy),2)+std::pow(((z_a-z_c)/ucz),2));
     
     
                 if (((d1 <= nnn_distance+0.01) && (d1 >= nnn_distance-0.01))&& ((d2 <= nnn_distance+0.01) && (d2 >= nnn_distance-0.01)) &&((d3 <= nnn_distance+0.01) && (d3 >= nnn_distance-0.01)) ){
                                          //get four spin exchange constant from material i to material j
                                          //add j k l to arrays for atom i .
                                          four_spin_neighbour_list_array_i.push_back(i);
                                          four_spin_neighbour_list_array_j.push_back(first_neigh[a]);
                                          four_spin_neighbour_list_array_k.push_back(first_neigh[b]);
                                          four_spin_neighbour_list_array_l.push_back(first_neigh[c]);

                                          four_spin_neighbour_list_array_i.push_back(first_neigh[a]);
                                          four_spin_neighbour_list_array_j.push_back(first_neigh[b]);
                                          four_spin_neighbour_list_array_k.push_back(first_neigh[c]);
                                          four_spin_neighbour_list_array_l.push_back(i);

                                          four_spin_neighbour_list_array_i.push_back(first_neigh[b]);
                                          four_spin_neighbour_list_array_j.push_back(first_neigh[c]);
                                          four_spin_neighbour_list_array_k.push_back(i);
                                          four_spin_neighbour_list_array_l.push_back(first_neigh[a]);
                                          
                                          four_spin_neighbour_list_array_i.push_back(first_neigh[c]);
                                          four_spin_neighbour_list_array_j.push_back(i);
                                          four_spin_neighbour_list_array_k.push_back(first_neigh[a]);
                                          four_spin_neighbour_list_array_l.push_back(first_neigh[b]);
                                          

                                    //    four times for all the permutations
                                           double fs_value= exchange::internal::mp[atoms::type_array[i]].fs[atoms::type_array[first_neigh[a]]];
                                          four_spin_exchange_list.push_back(fs_value*imus);
                                          four_spin_exchange_list.push_back(fs_value*imus);
                                          four_spin_exchange_list.push_back(fs_value*imus);
                                          four_spin_exchange_list.push_back(fs_value*imus);

                                          //4 interactions due to the permutations       
                                          n_interactions=n_interactions+1;                                    

                                          ofile << n_interactions << "\t" << i << "\t" <<  a << "\t" << b << "\t" << c <<"\t"<<-0.23e-21*imus<<std::endl;
                                          four_spin_neighbour_list_end_index[i] = n_interactions;
                                       }
                                    }
                                 }
                            }
                         }
                   
    ofile.close();
    std::cout<<"Four-spin quartets have been initialised"<<std::endl;

}
}
}
