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
   // double get_A(int material_index){
   //    if(internal::mp.empty()) return 0.0;
   //    if(material_index >= internal::mp.size()) return internal::mp[0].A.get(); // Fallback
   //    return internal::mp[material_index].A.get();
   // }

   // double get_Gamma(int material_index){
   //    if(internal::mp.empty()) return 0.0;
   //    if(material_index >= internal::mp.size()) return internal::mp[0].Gamma.get();
   //    return internal::mp[material_index].Gamma.get();
   // }

   // double get_omega0(int material_index){
   //    if(internal::mp.empty()) return 0.0;
   //    if(material_index >= internal::mp.size()) return internal::mp[0].omega0.get();
   //    return internal::mp[material_index].omega0.get();
   // }

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside quantum module
      //------------------------------------------------------------------------

      bool enabled = false; // bool to enable module

      std::vector<internal::mp_t> mp; // array of material properties

      noise_t noise_type = internal::quantum_zero; // Default to Quantum noise
      int window_size = 0; // Default window size (initially set to zero, later initialised to simulation length)
      int M_decimation = 1; // No decimation (interpolation in German -> English translation)
     

      std::vector<double> material_A_array; // compact array of Lorentzian amplitude for each material
      std::vector<double> material_gamma_array; // compact array of Lorentzian gamma for each material
      std::vector<double> material_omega0_array; // compact array of Lorentzian omega0 for each material

      std::vector<double> coarse_noise_field; // final generated noise (1D list of noise for all atoms and dimensions)

      std::vector<std::size_t> atom_idx_x; // indexing for accessing noise for each atom and spatial dimension in coarse array
      std::vector<std::size_t> atom_idx_y;
      std::vector<std::size_t> atom_idx_z;

      int noise_index = 0; // index for accessing noise arrays (at each fine time step)

      // Integration arrays for auxiliary variables
      std::vector<double> x_w_array;
      std::vector<double> y_w_array;
      std::vector<double> z_w_array;
      std::vector<double> x_v_array;
      std::vector<double> y_v_array;
      std::vector<double> z_v_array;

      // Integration arrays for RK4 integration
      std::vector<std::vector<double>> k1_storage;
      std::vector<std::vector<double>> k2_storage;
      std::vector<std::vector<double>> k3_storage;
      std::vector<std::vector<double>> k4_storage;
      std::vector<std::vector<double>> y_pred_storage;
      std::vector<std::vector<double>> y_in_storage;

      bool LLG_set = false; // delete this, use quantum::initialised

   } // end of internal namespace

} // end of quantum namespace
