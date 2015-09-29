#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

namespace program{


void setting_process(){

std::cout<< sim::H_applied << "\t" << sim::H_vec[0] << "\t" << sim::H_vec[1] << "\t" << sim::H_vec[2] <<std::endl;


double Configurations[8][3] = {{0,0,-1}, {-0.94,0,0.3}, {0.44, -0.845, 0.3},{0.44, 0.845, 0.3}, {0,0,1},{0.94,0,-0.3},{-0.44, 0.845, -0.3},{-0.44, -0.845, -0.3}};
double Possible_Arrays[8][4] = {{0,1,2,3},{1,0,3,2},{2,3,0,1},{3,2,1,0},{4,5,6,7},{5,4,7,6},{6,7,4,5},{7,6,5,4}};

		std::vector<int> No_in_Sublattice(4);

	// check calling of routine if error checking is activated
	if(err::check==true) std::cout << "program::setting has been called" << std::endl;

	for (int i = 0; i < atoms::num_atoms; i ++){
	//	std::cout<< i << "\t" << atoms::grain_array[i]<<"\t" << atoms::type_array[i] <<std::endl;
		for (int j = atoms::neighbour_list_start_index[i]; j < atoms::neighbour_list_end_index[i]; j ++){
				//std::cout<< i << "\t" << j << "\t" <<atoms::neighbour_list_array[j] <<"\t" << atoms::neighbour_interaction_type_array[j] << "\t" << atoms::type_array[atoms::neighbour_list_array[j]] << std::endl;
				if ((atoms::type_array[i] == 4) & (atoms::type_array[atoms::neighbour_list_array[j]] != 4)){
					No_in_Sublattice[atoms::type_array[atoms::neighbour_list_array[j]]] ++;
				}

			}
		}
double result,angle;
double min_angle = 1000;
double Direction_Closest_to_Field;
for (int i = 0; i < 8; i ++){

	//std::cout<< Configurations[i][0]*sim::H_vec[0] << "\t" << Configurations[i][1]*sim::H_vec[1] << "\t" << Configurations[i][2]*sim::H_vec[2] <<std::endl;
	result = Configurations[i][0]*sim::H_vec[0] + Configurations[i][1]*sim::H_vec[1] + Configurations[i][2]*sim::H_vec[2];
	angle = acos(result);
//	std::cout<<"angle:"<<angle << "\t" << result << "\t" <<std::endl;
	if ( angle< min_angle){
		Direction_Closest_to_Field = i;
		min_angle = angle;
	}
//	std::cout<<i << "\t" <<value<<"\t" <<min_angle <<std::endl;
}

std::cout<< No_in_Sublattice[0] << "\t" <<  No_in_Sublattice[1] << "\t" <<  No_in_Sublattice[2] << "\t" <<  No_in_Sublattice[3] << "\t" << std::endl;

int Max_atoms = 0;
int Largest_Sublattice;

for( int i = 0; i <4; i ++){
	if (No_in_Sublattice[i] > Max_atoms)
	{
		Largest_Sublattice = i;
		Max_atoms = No_in_Sublattice[i];
	}
}
int Chosen_array;
for (int i = 0; i <8; i ++){
//	std::cout<< Direction_Closest_to_Field<<"\t" <<Largest_Sublattice <<std::endl;
if (Possible_Arrays[i][Largest_Sublattice] == Direction_Closest_to_Field){
	Chosen_array = i;
	break;
}
//std::cout<< i <<std::endl;

}
//std::cout<< Chosen_array << std::endl;
std::cout<< Possible_Arrays[Chosen_array][0] <<"\t" << Possible_Arrays[Chosen_array][1] << "\t" << Possible_Arrays[Chosen_array][2] <<"\t" << Possible_Arrays[Chosen_array][3] <<std::endl;


for (int i = 0; i <atoms::num_atoms; i++){

		int Array = Possible_Arrays[Chosen_array][atoms::type_array[i]];
		atoms::x_spin_array[i] = Configurations[Array][0];
		atoms::y_spin_array[i] = Configurations[Array][1];
		atoms::z_spin_array[i] = Configurations[Array][2];

}

	// Calculate magnetisation statistics
	stats::mag_m();

	// Output data
	vout::data();
}


}//end of namespace program
