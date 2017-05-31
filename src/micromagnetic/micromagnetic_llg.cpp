
// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"

#include "random.hpp"
#include "errors.hpp"
#include "atoms.hpp"
#include "cells.hpp"
#include "sim.hpp"
#include "vmpi.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>

namespace mm = micromagnetic::internal;



namespace micromagnetic_arrays_llg{

	// Local arrays for LLG integration
	std::vector <double> x_euler_array;
	std::vector <double> y_euler_array;
	std::vector <double> z_euler_array;

	std::vector <double> x_array;
	std::vector <double> y_array;
	std::vector <double> z_array;

	std::vector <double> x_heun_array;
	std::vector <double> y_heun_array;
	std::vector <double> z_heun_array;

	std::vector <double> x_spin_storage_array;
	std::vector <double> y_spin_storage_array;
	std::vector <double> z_spin_storage_array;

	std::vector <double> x_initial_spin_array;
	std::vector <double> y_initial_spin_array;
	std::vector <double> z_initial_spin_array;

	bool LLG_set=false; ///< Flag to define state of LLG arrays (initialised/uninitialised)

}


namespace micromagnetic{

  int micromagnetic_init_llg(int num_cells){

  	// check calling of routine if error checking is activated
  	if(err::check==true) std::cout << "LLB_init has been called" << std::endl;

  	using namespace micromagnetic_arrays_llg;

	x_spin_storage_array.resize(num_cells,0.0);
  	y_spin_storage_array.resize(num_cells,0.0);
  	z_spin_storage_array.resize(num_cells,0.0);

	x_array.resize(num_cells,0.0);
  	y_array.resize(num_cells,0.0);
  	z_array.resize(num_cells,0.0);

  	x_initial_spin_array.resize(num_cells,0.0);
  	y_initial_spin_array.resize(num_cells,0.0);
  	z_initial_spin_array.resize(num_cells,0.0);

  	x_euler_array.resize(num_cells,0.0);
  	y_euler_array.resize(num_cells,0.0);
  	z_euler_array.resize(num_cells,0.0);

  	x_heun_array.resize(num_cells,0.0);
  	y_heun_array.resize(num_cells,0.0);
  	z_heun_array.resize(num_cells,0.0);

  	LLG_set=true;

   return EXIT_SUCCESS;
  }


int LLG( std::vector <int> local_cell_array,
							int num_steps,
                     int num_cells,
							int num_local_cells,
                     double temperature,
                     std::vector<double>& x_mag_array,
                     std::vector<double>& y_mag_array,
                     std::vector<double>& z_mag_array,
                     double Hx,
                     double Hy,
                     double Hz,
                     double H,
                     double dt,
                     std::vector <double> volume_array
                  ){

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "micromagnetic::LLG_Heun has been called" << std::endl;}

  using namespace micromagnetic_arrays_llg;

	// Check for initialisation of LLG integration arrays
	if(LLG_set== false) micromagnetic::micromagnetic_init_llg(num_cells);
	// Local variables for system integration

  //calculte chi(T).
	mm::one_o_chi_para =  mm::calculate_chi_para(num_local_cells,local_cell_array, num_cells, temperature);
	mm::one_o_chi_perp =  mm::calculate_chi_perp(num_local_cells,local_cell_array, num_cells, temperature);


  //The external fields equal the length of the field times the applied field vector.
  //This is saved to an array.
  mm::ext_field[0] = H*Hx;
  mm::ext_field[1] = H*Hy;
  mm::ext_field[2] = H*Hz;

  //save this new m as the initial value, so it can be saved and used in the final equation.
	for (int lc = 0; lc < num_local_cells; lc++){
		int cell = local_cell_array[lc];
		x_array[cell] = x_mag_array[cell];
		y_array[cell] = y_mag_array[cell];
		z_array[cell] = z_mag_array[cell];
    	double m_squared = sqrt(x_array[cell]*x_array[cell] + y_array[cell]*y_array[cell] + z_array[cell]*z_array[cell]);
    	x_array[cell] = x_array[cell]/m_squared;
    	y_array[cell] = y_array[cell]/m_squared;
    	z_array[cell] = z_array[cell]/m_squared;
		x_initial_spin_array[cell] = x_array[cell];
		y_initial_spin_array[cell] = y_array[cell];
		z_initial_spin_array[cell] = z_array[cell];
	}


  std::vector<double> m(3,0.0);
  std::vector<double> spin_field(3,0.0);


  for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){
		int cell = list_of_micromagnetic_cells[lc];
		m[0] = x_array[cell];
		m[1] = y_array[cell];
		m[2] = z_array[cell];
    //if (cell == 0) std::cout << m[0] << '\t'<< m[1] << '\t'<< m[2] << '\t' << field[0] <<'\t' << field[1] <<'\t' << field[2] << std::endl;
    spin_field = mm::calculate_llg_fields(m, temperature, num_cells, cell, x_array,y_array,z_array);

    const double one_oneplusalpha_sq = 1/(1+mm::alpha[cell]*mm::alpha[cell]); // material specific alpha and gamma
    const double alpha_oneplusalpha_sq = mm::alpha[cell]/(1+mm::alpha[cell]*mm::alpha[cell]);

    const double S[3] = {m[0],m[1],m[2]};
    const double H[3] = {-spin_field[0], -spin_field[1], -spin_field[2]};

    double xyz[3] = {0.0,0.0,0.0};
    xyz[0]=(one_oneplusalpha_sq)*(S[1]*H[2]-S[2]*H[1]) + (alpha_oneplusalpha_sq)*(S[1]*(S[0]*H[1]-S[1]*H[0])-S[2]*(S[2]*H[0]-S[0]*H[2]));
    xyz[1]=(one_oneplusalpha_sq)*(S[2]*H[0]-S[0]*H[2]) + (alpha_oneplusalpha_sq)*(S[2]*(S[1]*H[2]-S[2]*H[1])-S[0]*(S[0]*H[1]-S[1]*H[0]));
    xyz[2]=(one_oneplusalpha_sq)*(S[0]*H[1]-S[1]*H[0]) + (alpha_oneplusalpha_sq)*(S[0]*(S[2]*H[0]-S[0]*H[2])-S[1]*(S[1]*H[2]-S[2]*H[1]));

    // Store dS in euler array
    x_euler_array[cell] = xyz[0];
    y_euler_array[cell] = xyz[1];
    z_euler_array[cell] = xyz[2];

  }

  double S_new[3] = {0.0,0.0,0.0};
  double mod_S;

  for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){
		int cell = list_of_micromagnetic_cells[lc];

    S_new[0] = x_array[cell] + x_euler_array[cell]*dt;
    S_new[1] = y_array[cell] + y_euler_array[cell]*dt;
    S_new[2] = z_array[cell] + z_euler_array[cell]*dt;

   // Normalise Spin Length
   mod_S = 1.0/sqrt(S_new[0]*S_new[0] + S_new[1]*S_new[1] + S_new[2]*S_new[2]);

   S_new[0]=S_new[0]*mod_S;
   S_new[1]=S_new[1]*mod_S;
   S_new[2]=S_new[2]*mod_S;

   x_spin_storage_array[cell] = S_new[0];
   y_spin_storage_array[cell] = S_new[1];
   z_spin_storage_array[cell] = S_new[2];
//       if (cell == 0)std::cout << "store" <<x_spin_storage_array[cell] << std::endl;
 }

 for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){
  int cell = list_of_micromagnetic_cells[lc];
  m[0] = x_spin_storage_array[cell];
  m[1] = y_spin_storage_array[cell];
  m[2] = z_spin_storage_array[cell];

  spin_field = mm::calculate_llg_fields(m, temperature, num_cells, cell, x_array,y_array,z_array);
//  if (cell == 0)std::cout << m[0] << '\t' << field[0] << "\t" << field[1] << '\t' << field[2] << std::endl;

  const double S[3] = {m[0],m[1],m[2]};
  const double H[3] = {-spin_field[0], -spin_field[1], -spin_field[2]};


  const double one_oneplusalpha_sq = 1/(1+mm::alpha[cell]*mm::alpha[cell]); // material specific alpha and gamma
  const double alpha_oneplusalpha_sq = mm::alpha[cell]/(1+mm::alpha[cell]*mm::alpha[cell]);


  double xyz[3] = {0.0,0.0,0.0};
  xyz[0]=(one_oneplusalpha_sq)*(S[1]*H[2]-S[2]*H[1]) + (alpha_oneplusalpha_sq)*(S[1]*(S[0]*H[1]-S[1]*H[0])-S[2]*(S[2]*H[0]-S[0]*H[2]));
  xyz[1]=(one_oneplusalpha_sq)*(S[2]*H[0]-S[0]*H[2]) + (alpha_oneplusalpha_sq)*(S[2]*(S[1]*H[2]-S[2]*H[1])-S[0]*(S[0]*H[1]-S[1]*H[0]));
  xyz[2]=(one_oneplusalpha_sq)*(S[0]*H[1]-S[1]*H[0]) + (alpha_oneplusalpha_sq)*(S[0]*(S[2]*H[0]-S[0]*H[2])-S[1]*(S[1]*H[2]-S[2]*H[1]));

  // Store dS in euler array
   x_heun_array[cell] = xyz[0];
   y_heun_array[cell] = xyz[1];
   z_heun_array[cell] = xyz[2];
 }
 for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){
  int cell = list_of_micromagnetic_cells[lc];

  S_new[0]= x_initial_spin_array[cell]+(dt/2)*(x_euler_array[cell]+x_heun_array[cell]);
  S_new[1]= y_initial_spin_array[cell]+(dt/2)*(y_euler_array[cell]+y_heun_array[cell]);
  S_new[2]= z_initial_spin_array[cell]+(dt/2)*(z_euler_array[cell]+z_heun_array[cell]);

  // Normalise Spin Length
  mod_S = 1.0/sqrt(S_new[0]*S_new[0] + S_new[1]*S_new[1] + S_new[2]*S_new[2]);

  x_array[cell]=S_new[0]*mod_S;
  y_array[cell]=S_new[1]*mod_S;
  z_array[cell]=S_new[2]*mod_S;

  cells::mag_array_x[cell] = x_array[cell]*mm::ms[cell];
  cells::mag_array_y[cell] = y_array[cell]*mm::ms[cell];
  cells::mag_array_z[cell] = z_array[cell]*mm::ms[cell];

  }
  for(int atom_list=0;atom_list<number_of_none_atomistic_atoms;atom_list++){
		 int atom = list_of_none_atomistic_atoms[atom_list];
		 int cell = cells::atom_cell_id_array[atom];
		 atoms::x_spin_array[atom] = x_array[cell]*mm::m_e[cell];
		 atoms::y_spin_array[atom] = y_array[cell]*mm::m_e[cell];
		 atoms::z_spin_array[atom] = z_array[cell]*mm::m_e[cell];
     atoms::m_spin_array[atom] = mm::m_e[cell];

	}
  return 0;

}
}
