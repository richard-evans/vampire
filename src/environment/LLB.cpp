
// Vampire headers
#include "environment.hpp"

// micromagnetic module headers
#include "internal.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "cells.hpp"


#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>

namespace env = environment::internal;




namespace LLB_arrays{

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


namespace environment{

  int LLB_init(int num_cells){

		////std::cout << "called" <<std::endl;
  	using namespace LLB_arrays;
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

int LLB(double temperature,
        double Hx,
        double Hy,
        double Hz,
        double H,
        double dt){


  using namespace LLB_arrays;


  int num_cells = env::num_cells;

	// Check for initialisation of LLG integration arrays
	if(LLG_set== false) environment::LLB_init(num_cells);
	// Local variables for system integration

	//save this new m as the initial value, so it can be saved and used in the final equation.
	for (int i = 0; i < env::num_env_cells; i++){
		int cell = env::none_atomistic_cells[i];
		x_array[cell] = env::x_mag_array[cell]/env::Ms;
		y_array[cell] = env::y_mag_array[cell]/env::Ms;
		z_array[cell] = env::z_mag_array[cell]/env::Ms;
		x_initial_spin_array[cell] = x_array[cell];
		y_initial_spin_array[cell] = y_array[cell];
		z_initial_spin_array[cell] = z_array[cell];
	}
//	std::cout << 	x_array[0] << '\t' << y_array[0] << '\t' << z_array[0] << std::endl;
   //if (sim::time % demag_update_rate == 0) env::calculate_demag_fields();
	const double kB = 1.3806503e-23;
	std::vector<double> m(3,0.0);
	std::vector<double> spin_field(3,0.0);

  env::ext_field[0] = H*Hx;
  env::ext_field[1] = H*Hy;
	env::ext_field[2] = H*Hz;

  //calculte chi(T).
	env::one_o_chi_para =  env::calculate_chi_para(temperature);
	env::one_o_chi_perp =  env::calculate_chi_perp(temperature);

  //6 arrays of gaussian random numbers to store the stochastic noise terms for x,y,z parallel and perperdicular
	std::vector <double> GW1x(num_cells,0.0);
	std::vector <double> GW1y(num_cells,0.0);
	std::vector <double> GW1z(num_cells,0.0);
	std::vector <double> GW2x(num_cells,0.0);
	std::vector <double> GW2y(num_cells,0.0);
	std::vector <double> GW2z(num_cells,0.0);

//std::cout << "be" <<std::endl;

	//fill the noise terms
	generate (GW1x.begin(),GW1x.end(), mtrandom::gaussian);
	generate (GW1y.begin(),GW1y.end(), mtrandom::gaussian);
	generate (GW1z.begin(),GW1z.end(), mtrandom::gaussian);
	generate (GW2x.begin(),GW2x.end(), mtrandom::gaussian);
	generate (GW2y.begin(),GW2y.end(), mtrandom::gaussian);
	generate (GW2z.begin(),GW2z.end(), mtrandom::gaussian);

	for (int i = 0; i < env::num_env_cells; i++){
		int cell = env::none_atomistic_cells[i];

    m[0] = x_array[cell];
		m[1] = y_array[cell];
		m[2] = z_array[cell];
		   if (cell == 0 && sim::time % demag_update_rate == 0) env::calculate_demag_fields();
    spin_field = env::calculate_llb_fields(m, temperature, cell, x_array,y_array,z_array);

    double sigma_para = sqrt(2*kB*temperature*env::alpha_para/(env::Ms*dt)); //why 1e-27
    double sigma_perp = sqrt(2*kB*temperature*(env::alpha_perp-env::alpha_para)/(dt*env::Ms*env::alpha_perp*env::alpha_perp));
    const double H[3] = {spin_field[0], spin_field[1], spin_field[2]};

    //saves the noise terms to an array
    const double GW2t[3] = {GW2x[cell],GW2y[cell],GW2z[cell]};
    const double one_o_m_squared = 1.0/(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
    const double SdotH = m[0]*H[0] + m[1]*H[1] + m[2]*H[2];
    double xyz[3] = {0.0,0.0,0.0};
		//calculates the LLB equation
		xyz[0]=  - (m[1]*H[2]-m[2]*H[1])
						 + env::alpha_para*m[0]*SdotH*one_o_m_squared
						 - env::alpha_perp*(m[1]*(m[0]*H[1]-m[1]*H[0])-m[2]*(m[2]*H[0]-m[0]*H[2]))*one_o_m_squared
						 + GW1x[cell]*sigma_para
						 - env::alpha_perp*(m[1]*(m[0]*GW2t[1]-m[1]*GW2t[0])-m[2]*(m[2]*GW2t[0]-m[0]*GW2t[2]))*one_o_m_squared*sigma_perp;

		xyz[1]=  - (m[2]*H[0]-m[0]*H[2])
						 + env::alpha_para*m[1]*SdotH*one_o_m_squared
						 - env::alpha_perp*(m[2]*(m[1]*H[2]-m[2]*H[1])-m[0]*(m[0]*H[1]-m[1]*H[0]))*one_o_m_squared
						 + GW1y[cell]*sigma_para
						 - env::alpha_perp*(m[2]*(m[1]*GW2t[2]-m[2]*GW2t[1])-m[0]*(m[0]*GW2t[1]-m[1]*GW2t[0]))*one_o_m_squared*sigma_perp;

		xyz[2]=	 - (m[0]*H[1]-m[1]*H[0])
						 + env::alpha_para*m[2]*SdotH*one_o_m_squared
						 - env::alpha_perp*(m[0]*(m[2]*H[0]-m[0]*H[2])-m[1]*(m[1]*H[2]-m[2]*H[1]))*one_o_m_squared
						 + GW1z[cell]*sigma_para
						 - env::alpha_perp*(m[0]*(m[2]*GW2t[0]-m[0]*GW2t[2])-m[1]*(m[1]*GW2t[2]-m[2]*GW2t[1]))*one_o_m_squared*sigma_perp;

		x_euler_array[cell] = xyz[0];
		y_euler_array[cell] = xyz[1];
		z_euler_array[cell] = xyz[2];

	}

  for (int cell = 0; cell < num_cells; cell++){
    x_spin_storage_array[cell] = x_array[cell] + x_euler_array[cell]*dt;
    y_spin_storage_array[cell] = y_array[cell] + y_euler_array[cell]*dt;
    z_spin_storage_array[cell] = z_array[cell] + z_euler_array[cell]*dt;
  }


	for (int i = 0; i < env::num_env_cells; i++){
		int cell = env::none_atomistic_cells[i];
    m[0] = x_array[cell];
		m[1] = y_array[cell];
		m[2] = z_array[cell];
    spin_field = env::calculate_llb_fields(m, temperature, cell, x_array,y_array,z_array);

    double sigma_para = sqrt(2*kB*temperature*env::alpha_para/(env::Ms*dt)); //why 1e-27
    double sigma_perp = sqrt(2*kB*temperature*(env::alpha_perp-env::alpha_para)/(dt*env::Ms*env::alpha_perp*env::alpha_perp));
    const double H[3] = {spin_field[0], spin_field[1], spin_field[2]};

    //saves the noise terms to an array
    const double GW2t[3] = {GW2x[cell],GW2y[cell],GW2z[cell]};
    const double one_o_m_squared = 1.0/(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
    const double SdotH = m[0]*H[0] + m[1]*H[1] + m[2]*H[2];
    double xyz[3] = {0.0,0.0,0.0};
		//calculates the LLB equation
		xyz[0]=  - (m[1]*H[2]-m[2]*H[1])
						 + env::alpha_para*m[0]*SdotH*one_o_m_squared
						 - env::alpha_perp*(m[1]*(m[0]*H[1]-m[1]*H[0])-m[2]*(m[2]*H[0]-m[0]*H[2]))*one_o_m_squared
						 + GW1x[cell]*sigma_para
						 - env::alpha_perp*(m[1]*(m[0]*GW2t[1]-m[1]*GW2t[0])-m[2]*(m[2]*GW2t[0]-m[0]*GW2t[2]))*one_o_m_squared*sigma_perp;

		xyz[1]=  - (m[2]*H[0]-m[0]*H[2])
						 + env::alpha_para*m[1]*SdotH*one_o_m_squared
						 - env::alpha_perp*(m[2]*(m[1]*H[2]-m[2]*H[1])-m[0]*(m[0]*H[1]-m[1]*H[0]))*one_o_m_squared
						 + GW1y[cell]*sigma_para
						 - env::alpha_perp*(m[2]*(m[1]*GW2t[2]-m[2]*GW2t[1])-m[0]*(m[0]*GW2t[1]-m[1]*GW2t[0]))*one_o_m_squared*sigma_perp;

		xyz[2]=	 - (m[0]*H[1]-m[1]*H[0])
						 + env::alpha_para*m[2]*SdotH*one_o_m_squared
						 - env::alpha_perp*(m[0]*(m[2]*H[0]-m[0]*H[2])-m[1]*(m[1]*H[2]-m[2]*H[1]))*one_o_m_squared
						 + GW1z[cell]*sigma_para
						 - env::alpha_perp*(m[0]*(m[2]*GW2t[0]-m[0]*GW2t[2])-m[1]*(m[1]*GW2t[2]-m[2]*GW2t[1]))*one_o_m_squared*sigma_perp;

		x_heun_array[cell] = xyz[0];
		y_heun_array[cell] = xyz[1];
		z_heun_array[cell] = xyz[2];

	}

	for (int i = 0; i < env::num_env_cells; i++){
		int cell = env::none_atomistic_cells[i];
    x_array[cell] = x_initial_spin_array[cell] + 0.5*dt*(x_euler_array[cell] + x_heun_array[cell]);
    y_array[cell] = y_initial_spin_array[cell] + 0.5*dt*(y_euler_array[cell] + y_heun_array[cell]);
    z_array[cell] = z_initial_spin_array[cell] + 0.5*dt*(z_euler_array[cell] + z_heun_array[cell]);


    env::x_mag_array[cell] = x_array[cell]*env::Ms;
    env::y_mag_array[cell] = y_array[cell]*env::Ms;
    env::z_mag_array[cell] = z_array[cell]*env::Ms;
  }




if (sim::time %1000 == 0) 	int a = env::output();


	return 0;

}


}
