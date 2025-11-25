//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Fried-Conrad Weber 2025. All rights reserved.
//
//   Email: fried-conrad.weber@uni-potsdam.de
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "quantum.hpp"

// quantum module headers
#include "internal.hpp"

namespace quantum{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   //---------------------------------------------------------------------------
   // Material parameter accessors
   //---------------------------------------------------------------------------
   double get_A(int material_index){
       if(internal::mp.empty()) return 0.0;
       if(material_index >= internal::mp.size()) return internal::mp[0].A.get(); // Fallback
       return internal::mp[material_index].A.get();
   }

   double get_Gamma(int material_index){
       if(internal::mp.empty()) return 0.0;
       if(material_index >= internal::mp.size()) return internal::mp[0].Gamma.get();
       return internal::mp[material_index].Gamma.get();
   }

   double get_omega0(int material_index){
       if(internal::mp.empty()) return 0.0;
       if(material_index >= internal::mp.size()) return internal::mp[0].omega0.get();
       return internal::mp[material_index].omega0.get();
   }

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside quantum module
      //------------------------------------------------------------------------

      bool enabled = false; // bool to enable module

      std::vector<internal::mp_t> mp; // array of material properties

      int noise_type = 1; // Default to Quantum noise
      int window_size = 6000; // Default window size
      int M_decimation = 1;
      double noise_index = 0.0;

      std::vector<double> coarse_noise_field;
      std::vector<std::size_t> atom_idx_x;
      std::vector<std::size_t> atom_idx_y;
      std::vector<std::size_t> atom_idx_z;

      // Integration arrays
      std::vector<double> x_w_array;
      std::vector<double> y_w_array;
      std::vector<double> z_w_array;
      std::vector<double> x_v_array;
      std::vector<double> y_v_array;
      std::vector<double> z_v_array;

      std::vector<std::vector<double>> k1_storage;
      std::vector<std::vector<double>> k2_storage;
      std::vector<std::vector<double>> k3_storage;
      std::vector<std::vector<double>> k4_storage;
      std::vector<std::vector<double>> y_pred_storage;
      std::vector<std::vector<double>> y_in_storage;

      bool LLG_set = false;

   } // end of internal namespace

} // end of quantum namespace

