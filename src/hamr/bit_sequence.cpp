//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) Andrea Meo 2022. 
//
// All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <iostream>
// #include "math.h"

// Vampire headers
// #include "material.hpp"
// #include "vmpi.hpp"
#include "errors.hpp"
#include "hamr.hpp"
#include "vio.hpp"

// hamr headers
#include "internal.hpp"

namespace hamr{
   namespace internal{

      //-----------------------------------------------------------------------------
      // Function to generate single tone sequence
      //-----------------------------------------------------------------------------
		void create_singletone_vector(){

			if(!hamr::internal::create_singletone){
				std::cout << "Bit sequence provided by user" << std::endl;
				zlog << zTs() << "Bit sequence provided by user" << std::endl;
				return;
			}

			// output informative message
			zlog << zTs() << "Creating singletone bit sequence ..." << std::endl;

			hamr::internal::bit_sequence.clear();
			for(auto i=0; i<hamr::internal::num_bits; ++i){
				int bit = i%2;
				if(bit==0){ bit=-1;}
				hamr::internal::bit_sequence.push_back(bit);
			}

			return;
		}
      

      //-----------------------------------------------------------------------------
      // Function to check that user defined bit sequence is consistent with system
		// size and "number-of-bits" parameter
      //-----------------------------------------------------------------------------
		void check_sequence_length(){

			zlog << zTs() << "Checking length of bit sequence ..." << std::endl;

      	// Check that "number-of-bits" and size of bit sequence provided are consistent
      	if(hamr::internal::num_bits > hamr::internal::bit_sequence.size()){
      	   std::cout << "Warning: Requested number-of-bits "  << hamr::internal::num_bits 
      	            << " larger than size of the provided bit sequence=" << hamr::internal::bit_sequence.size() 
      	            << ". Adjusting to " << hamr::internal::bit_sequence.size() << std::endl;
      	   zlog << zTs() << "Warning: Requested number-of-bit "  << hamr::internal::num_bits 
      	               << " larger than size of the provided bit sequence=" << hamr::internal::bit_sequence.size() 
      	               << ". Adjusting to " << hamr::internal::bit_sequence.size() << std::endl;
      	   hamr::internal::num_bits = hamr::internal::bit_sequence.size();
      	}
      	else if(hamr::internal::num_bits < hamr::internal::bit_sequence.size()){ 
      	   std::cout << "Warning: number-of-bits "  << hamr::internal::num_bits 
      	            << " smaller than size of provided bit sequence=" << hamr::internal::bit_sequence.size() 
      	            << ". Trimming bit sequence." << std::endl;
      	   zlog << zTs() << "Warning: number-of-bits "  << hamr::internal::num_bits 
      	               << " smaller than size of provided bit sequence=" << hamr::internal::bit_sequence.size() 
      	               << ". Trimming bit sequence." << std::endl;
				hamr::internal::bit_sequence.resize(hamr::internal::num_bits);
      	}

      	// Check that number of bit requested is compatible with system size
      	if(hamr::internal::num_bits > hamr::internal::num_tracks * hamr::internal::bits_per_tack){
				const int num_bits_total = hamr::internal::num_tracks * hamr::internal::bits_per_tack;
				std::cout << "Warning: number-of-bits "  << hamr::internal::num_bits 
				         << " too big for system size. Reducing to " << num_bits_total << std::endl;
				zlog << zTs() << "Warning: number-of-bits "  << hamr::internal::num_bits 
				               << " too big for system size. Reducing to " << num_bits_total << std::endl;
				hamr::internal::num_bits = num_bits_total;
				hamr::internal::bit_sequence.resize(num_bits_total);
      	}

			return;
		}


   } // end of namespace internal
} // end of namespace hamr
