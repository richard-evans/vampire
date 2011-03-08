#include "atoms.hpp"
#include "errors.hpp"
#include "program.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include <iostream>

namespace program{
	
int bmark(){
  // check calling of routine if error checking is activated
  if(err::check==true){std::cout << "program::bmark has been called" << std::endl;}

  // Setup LLG arrays
  sim::initialise();

	// Initialise spins to random state
	for(int atom =0;atom<atoms::num_atoms;atom++){
		atoms::x_spin_array[atom]=0.0;
		atoms::y_spin_array[atom]=0.0;
		atoms::z_spin_array[atom]=1.0;
  }

  sim::temperature=300.0;

  	sim::integrator=0;
  
  // Simulate system
  for(sim::time=0;sim::time<sim::total_time;sim::time+=sim::partial_time){

  // Calculate LLG
  //sim::LLG(1);
  	sim::integrate(sim::partial_time);


      // Calculate mag_m, mag
  //if(sim::time%sim::partial_time==0){
      stats::mag_m();
		//vout::pov_file();
  if(vmpi::my_rank==0){
    std::cout << sim::time << "\t" << stats::total_mag_m_norm;
    std::cout << "\t" << stats::total_mag_norm[0];
    std::cout << "\t" << stats::total_mag_norm[1];
    std::cout << "\t" << stats::total_mag_norm[2];
    std::cout << std::endl;
    //vmag << sim::temperature << "\t" << stats::total_mag_m_norm << std::endl;
  }

  }

return EXIT_SUCCESS;
}

}//end of namespace program

