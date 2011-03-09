#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include <iostream>

namespace program{

/// @brief Function to calculate the hysteresis loop
///
/// @callgraph
/// @callergraph
///
/// @details Consists of a sequence of sub-calculations of fixed temperature. The system is initialised 
/// ordered. After initialisation a whole hysteresis loop of the system and coercivity are calculated.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Weijia Fan, wf507@york.ac.uk
/// @version 1.0
/// @date    27/01/2010
///
/// @return EXIT_SUCCESS
/// 
/// @internal
///	Created:		27/01/2010
///	Revision:	  ---
///=====================================================================================
///
int hysteresis(){
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "hysteresis has been called" << std::endl;}
	
	// Declare function prototype
	//int calculate_applied_fields(const int,const int);
	// set loop temperature
	
	// Initialise spins to +z
	for(int atom =0;atom<atoms::num_atoms;atom++){
		atoms::x_spin_array[atom]=0.0;
		atoms::y_spin_array[atom]=0.0;
		atoms::z_spin_array[atom]=1.0;
		//atoms::z_spin_array[atom]=1.0-2.0*double(vmpi::my_rank%2);
	}
	
	 // Setup Hmax and J=Hinc
	 int Hmax = 6000; // mT
	 int Hinc = 10; // mT	 
	 
	 // setup mag_new(H_new) and mag_old(H_old) for coercive field
	 double mag_new = 0.0;
	 double mag_old = 0.0;
	 double H_old = 0.0;
	 double H_new = 0.0;
	 double H_c_left= 0.0;
	 double H_c_right = 0.0;
	 double H_coercive = 0.0;
	 //std::cout << std::cout.flags() << std::endl;
        std::cout.unsetf(std::ios_base::floatfield);
        //std::cout << std::cout.flags() << std::endl;
	 //std::cout << mp::material[0].mu_s_SI << "\t" << mp::material[0].Ku1_SI << "\t" << 2.0*mp::material[0].Ku1_SI/mp::material[0].mu_s_SI << std::endl;
	 std::cout << "Estimated Coercivity:" << 2.0*mp::material[0].Ku1_SI/mp::material[0].mu_s_SI << std::endl;
	 vout::pov_file();
	 // parity loop
	 for(int parity=-1;parity<2;parity+=2){
	  // Set up loop variables
	  for (int H = -Hmax; H<= Hmax;H+=Hinc){

	    sim::H_applied=double(parity)*double(H)*0.001;	// Tesla
	  
	  // time loop
	    for(sim::time=0;sim::time<sim::loop_time;sim::time+=sim::partial_time){
	      sim::LLG(sim::partial_time);
	      stats::mag_m();
	      //if(vmpi::my_rank==0){
	      // 	std::cout << sim::time<< "\t" << stats::total_mag_m_norm;
	      //	std::cout << "\t" << stats::total_mag_norm[0];
	      //	std::cout << "\t" << stats::total_mag_norm[1];
	      //	std::cout << "\t" << stats::total_mag_norm[2];
	      //	std::cout << std::endl;
	      //}
	    } 
	    
		// output pov_file after each field point
	    //vout::pov_file();
		
	    std::cout << vmpi::my_rank;
	    std::cout << "\t" << stats::total_mag_m_norm;
	    std::cout << "\t" << stats::total_mag_norm[0];
	    std::cout << "\t" << stats::total_mag_norm[1];
	    std::cout << "\t" << stats::total_mag_norm[2];
	    std::cout << std::endl;
	    
	     if(vmpi::my_rank==0){
		mag_new = stats::total_mag_norm[2];
		H_new = sim::H_applied;
		if ((mag_new*mag_old < 0) && (mag_old > mag_new)){
		  // calculate coercive field
		  H_c_left=(mag_old*H_new-mag_new*H_old)/(mag_old-mag_new);
		  std::cout << "\t" << "the left coercivity is" << "\t" << H_c_left << "\t" << "Tesla" << std::endl;	
		}
		if ((mag_new*mag_old < 0) && (mag_old < mag_new)){
		  H_c_right=(mag_old*H_new-mag_new*H_old)/(mag_old-mag_new);
		  std::cout << "\t" << "the right coercivity is" << "\t" <<  H_c_right << "\t" << "Tesla" << std::endl;
		}		 		  
	      std::cout << sim::H_applied << "\t" << stats::total_mag_m_norm; // Tesla
	      std::cout << "\t" << stats::total_mag_norm[0];
	      std::cout << "\t" << stats::total_mag_norm[1];
	      std::cout << "\t" << stats::total_mag_norm[2];
	      std::cout << std::endl;
		// check current and before values
		//std::cout << "current mag_new is" << mag_new << "\t" << "before mag_old is" << mag_old<< std::endl;
		//std::cout << "current H_new is" << H_new << "\t" << "before H is" << H_old << std::endl;
	      vmag << sim::H_applied << "\t" << stats::total_mag_m_norm;
	      vmag << "\t" << stats::total_mag_norm[0];
	      vmag << "\t" << stats::total_mag_norm[1];
	      vmag << "\t" << stats::total_mag_norm[2];
	      vmag << std::endl;
	      

	      				
	      mag_old = mag_new;
	      H_old = H_new;
	    }
		 

	  //if ((H%100)==0){vout::pov_file();}		
	     // mag_new = stats::total_mag_norm[2];
	    
	     // if ((mag_old-0.8)/(mag_new-0.8)<0){
		//vout::pov_file();	      
	     // }
	     // if ((mag_old+0.8)/(mag_new+0.8)<0){
		//vout::pov_file();	      
	     // }
      	     // mag_old = mag_new;

	  }


	 }
	 if(vmpi::my_rank==0){
	  H_coercive = -H_c_left;   //0.5*(H_c_right-H_c_left);
	  std::cout << "coercive field of the system is" << "\t" << H_coercive << "\t" << "Tesla" << std::endl;
	  //vmag << "Hc+ is" << "\t" << H_c_right << "\tTesla" << "Hc- is" << "\t" << H_c_left << "\tTesla" << std::endl;
	  vmag << "The coercive field of the system is" << "\t" << H_coercive << "\t" << "Tesla" << std::endl;
	 }
	 
	return EXIT_SUCCESS;
  }

}//end of namespace program


