//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2011-2018. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
#ifndef STATS_H_
#define STATS_H_

// C++ include files
#include <vector>
#include <string>

namespace stats
//==========================================================
// Namespace statistics
//==========================================================
{
	extern int num_atoms;				//Number of atoms for statistic purposes

	// Member Functions
	extern double max_torque();

   /// Statistics energy types
   enum energy_t { total = 0, exchange = 1, anisotropy = 2, applied_field = 3, magnetostatic = 4};

   /// Statistics types
   enum stat_t { atotal=0, mean=1};

   //-------------------------------------------------
   // New statistics module functions and variables
   //-------------------------------------------------

   // Control functions
   void initialize(const int num_atoms,
                   const int num_materials,
                   const int num_grains,
                   const std::vector<double>& magnetic_moment_array,
                   const std::vector<int>& material_type_array,
                   const std::vector<int>& grain_array,
                   const std::vector<int>& height_category_array,
                   const std::vector<bool>& non_magnetic_materials_array);


	// Function to update statistics
	void update();

	// Function to reset average statistics counters
   void reset();

	// Statistics control flags (to be moved internally when long-awaited refactoring of vio is done)
	extern bool calculate_system_energy;
	extern bool calculate_grain_energy;
	extern bool calculate_material_energy;

	extern bool calculate_system_magnetization;
	extern bool calculate_grain_magnetization;
	extern bool calculate_material_magnetization;
	extern bool calculate_material_grain_magnetization;
	extern bool calculate_height_magnetization;
	extern bool calculate_material_height_magnetization;
	extern bool calculate_material_grain_height_magnetization;

	extern bool calculate_system_torque;
	extern bool calculate_grain_torque;
	extern bool calculate_material_torque;

	extern bool calculate_system_specific_heat;
	extern bool calculate_grain_specific_heat;
	extern bool calculate_material_specific_heat;

	extern bool calculate_material_standard_deviation;

	extern bool calculate_system_susceptibility;
	extern bool calculate_grain_susceptibility;
	extern bool calculate_material_susceptibility;

	extern bool calculate_system_binder_cumulant;
	extern bool calculate_material_binder_cumulant;

   extern bool calculate_system_spin_length;
   extern bool calculate_material_spin_length;
   extern bool calculate_height_spin_length;

	// forward declaration of friend classes
	class susceptibility_statistic_t;
	class specific_heat_statistic_t;
        class binder_cumulant_statistic_t;

	class standard_deviation_statistic_t;
   //----------------------------------
   // Energy class definition
   //----------------------------------
   class energy_statistic_t{

      friend class specific_heat_statistic_t;

   public:
      energy_statistic_t (std::string n):initialized(false){
        name = n;
      };
      bool is_initialized();
      void set_mask(const int in_mask_size, const std::vector<int> in_mask);
      void get_mask(std::vector<int>& out_mask, std::vector<double>& out_normalisation);
      void calculate(const std::vector<double>& sx, const std::vector<double>& sy, const std::vector<double>& sz,
                     const std::vector<double>& mm, const std::vector<int>& mat, const double temperature);

      void reset_averages();

      void set_total_energy(         std::vector<double>& new_energy, std::vector<double>& new_mean_energy);
      void set_exchange_energy(      std::vector<double>& new_energy, std::vector<double>& new_mean_energy);
      void set_anisotropy_energy(    std::vector<double>& new_energy, std::vector<double>& new_mean_energy);
      void set_applied_field_energy( std::vector<double>& new_energy, std::vector<double>& new_mean_energy);
      void set_magnetostatic_energy( std::vector<double>& new_energy, std::vector<double>& new_mean_energy);

      const std::vector<double>& get_total_energy();
      const std::vector<double>& get_exchange_energy();
      const std::vector<double>& get_anisotropy_energy();
      const std::vector<double>& get_applied_field_energy();
      const std::vector<double>& get_magnetostatic_energy();

      void update_mean_counter(long counter);

      std::string output_energy(enum energy_t energy_type, bool header);
      std::string output_mean_energy(enum energy_t energy_type, bool header);

   private:
      bool initialized;
      int num_atoms;
      int mask_size;
      double mean_counter;

      std::vector<int> mask;

      std::vector<double> total_energy;
      std::vector<double> exchange_energy;
      std::vector<double> anisotropy_energy;
      std::vector<double> applied_field_energy;
      std::vector<double> magnetostatic_energy;

      std::vector<double> mean_total_energy;
      std::vector<double> mean_exchange_energy;
      std::vector<double> mean_anisotropy_energy;
      std::vector<double> mean_applied_field_energy;
      std::vector<double> mean_magnetostatic_energy;

      std::vector<int> zero_list;
      std::vector<double> normalisation;

      std::string name;

   };

   //----------------------------------
   // Magnetization Class definition
   //----------------------------------
   class magnetization_statistic_t{

      friend class susceptibility_statistic_t;
      friend class standard_deviation_statistic_t;
      friend class spin_length_statistic_t;
      friend class binder_cumulant_statistic_t;
      public:
         magnetization_statistic_t (std::string n):initialized(false){
           name = n;
         };
         bool is_initialized();
         void set_mask(const int mask_size, std::vector<int> inmask, const std::vector<double>& mm);
         void get_mask(std::vector<int>& out_mask, std::vector<double>& out_saturation);
         void calculate_magnetization(const std::vector<double>& sx, const std::vector<double>& sy, const std::vector<double>& sz, const std::vector<double>& mm);
         void set_magnetization(std::vector<double>& magnetization, std::vector<double>& mean_magnetization, long counter);
         void reset_magnetization_averages();
         const std::vector<double>& get_magnetization();
         void save_checkpoint(std::ofstream& chkfile);
         void load_checkpoint(std::ifstream& chkfile, bool chk_continue);
         const std::vector<double>& get_checkpoint_parameters(double& sum_mx, double& sum_my, double& sum_mz, double& sum_count);
         std::string output_magnetization(bool header);
         std::string output_normalized_magnetization(bool header);
         std::string output_normalized_magnetization_length(bool header);
         std::string output_normalized_mean_magnetization(bool header);
         std::string output_normalized_mean_magnetization_length(bool header);
         std::string output_normalized_magnetization_dot_product(const std::vector<double>& vec,bool header);
         std::string output_mean_magnetization_length(bool header);
         std::string output_mean_magnetization(bool header);

      private:
         bool initialized;
         int num_atoms;
         int mask_size;
         double mean_counter;
         std::vector<int> mask;
         std::vector<double> magnetization;
         std::vector<double> mean_magnetization;
         std::vector<int> zero_list;
         std::vector<double> saturation;
         std::string name;

   };

	//----------------------------------
   // Torque class definition
   //----------------------------------
   class torque_statistic_t{

      public:
         torque_statistic_t (std::string n):initialized(false){
           name = n;
         };
         bool is_initialized();
         void set_mask(const int mask_size, std::vector<int> inmask, const std::vector<double>& mm);
         void get_mask(std::vector<int>& out_mask);
         void calculate_torque(const std::vector<double>& sx, const std::vector<double>& sy, const std::vector<double>& sz,
										 const std::vector<double>& bxs, const std::vector<double>& bys, const std::vector<double>& bzs,
										 const std::vector<double>& bxe, const std::vector<double>& bye, const std::vector<double>& bze,
										 const std::vector<double>& mm);
         void set_torque(std::vector<double>& torque, std::vector<double>& mean_torque, long counter);
         void reset_torque_averages();
         const std::vector<double>& get_torque();
         std::string output_torque(bool header);
			std::string output_mean_torque(bool header);

      private:
         bool initialized;
         int num_atoms;
         int mask_size;
         double mean_counter;
         std::vector<int> mask;
			std::vector<int> num_atoms_in_mask;
         std::vector<double> torque;
         std::vector<double> mean_torque;
         std::vector<int> zero_list;
         std::string name;

   };

	//----------------------------------
   // Specific Heat Class definition
   //----------------------------------
   class specific_heat_statistic_t{

      public:
			specific_heat_statistic_t (std::string n):initialized(false){
				name = n;
			};
			void initialize(energy_statistic_t& energy_statistic);
			void calculate(const std::vector<double>& energy);
			void save_checkpoint(std::ofstream& chkfile);
			void load_checkpoint(std::ifstream& chkfile, bool chk_continue);
			void reset_averages();
			std::string output_mean_specific_heat(const double temperature,bool header);


      private:
         bool initialized;
         int num_elements;
         double mean_counter;
         std::vector<double> mean_specific_heat;
         std::vector<double> mean_specific_heat_squared;
         std::vector<double> normalisation;

         std:: string name;

   };

   //----------------------------------
   // Susceptibility Class definition
   //----------------------------------
   class susceptibility_statistic_t{

      public:
         susceptibility_statistic_t (std::string n):initialized(false){
           name = n;
         };
			void initialize(magnetization_statistic_t& mag_stat);
			void calculate(const std::vector<double>& magnetization);
			void save_checkpoint(std::ofstream& chkfile);
			void load_checkpoint(std::ifstream& chkfile, bool chk_continue);
			void reset_averages();
			std::string output_mean_susceptibility(const double temperature,bool header);
         //std::string output_mean_absolute_susceptibility();

      private:
         bool initialized;
         int num_elements;
         double mean_counter;
         std::vector<double> mean_susceptibility;
         std::vector<double> mean_susceptibility_squared;
         std::vector<double> mean_absolute_susceptibility;
         std::vector<double> mean_absolute_susceptibility_squared;
         std::vector<double> saturation;
         std::string name;

   };

   //----------------------------------
   // Spin Length Class definition
   //----------------------------------
   class spin_length_statistic_t{

      public:
         spin_length_statistic_t (std::string n):initialized(false){
           name = n;
         };
         bool is_initialized();
         void set_mask(const int mask_size, std::vector<int> inmask);
         void get_mask(std::vector<int>& out_mask);
         void calculate_spin_length(const std::vector<double>& sx, const std::vector<double>& sy, const std::vector<double>& sz);
         void reset_averages();
         std::string output_mean_spin_length(bool header);

      private:
         bool initialized;
         int num_atoms;
         int mask_size;
         double mean_counter;
         std::vector<int> mask;
         std::vector<double> spin_length;
         std::vector<double> mean_spin_length;
         std::vector<int> zero_list;
         std::vector<double> normalisation;
         std::string name;

   };

   //----------------------------------
   // Standard Deviation of magnetisation in time Class definition
   //----------------------------------
   class standard_deviation_statistic_t{ // AJN

      public:
         standard_deviation_statistic_t (std::string n):initialized(false){
           name=n;
         };
         void initialize(magnetization_statistic_t& mag_stat);
         void update(const std::vector<double>& magnetization);
         void reset_averages();
         std::string output_standard_deviation(bool header);

      private:
         bool initialized;
         int num_elements;//number of elements in the system
         int idx;// index for looping through directions
         double mean_counter;// counts time steps in loop - factor out for normal time
         double res1; // residuals calculated at each time
         double res2;
         std::vector<double> residual_sq;// running squared residual for each direction
         std::vector<double> mean; // running mean for each direction
         std::string name;

   };
   //----------------------------------
   // Binder cumulant Class definition
   //----------------------------------
   class binder_cumulant_statistic_t{

      public:
         binder_cumulant_statistic_t (std::string n):initialized(false){
           name = n;
         };
         void initialize(magnetization_statistic_t& mag_stat);
         void calculate(const std::vector<double>& magnetization);
         void reset_averages();
         std::string output_binder_cumulant(bool header);

      private:
         bool initialized;
         int num_elements;
         double mean_counter;
         std::vector<double> binder_cumulant_squared;
         std::vector<double> binder_cumulant_fourth_power;
         std::string name;

   };

   //----------------------------------
	// Statistics class instantiations
   //----------------------------------
	extern energy_statistic_t system_energy;
	extern energy_statistic_t grain_energy;
	extern energy_statistic_t material_energy;

   extern magnetization_statistic_t system_magnetization;
	extern magnetization_statistic_t grain_magnetization;
	extern magnetization_statistic_t material_magnetization;
	extern magnetization_statistic_t material_grain_magnetization;
   extern magnetization_statistic_t height_magnetization;
   extern magnetization_statistic_t material_height_magnetization;
   extern magnetization_statistic_t material_grain_height_magnetization;

	extern torque_statistic_t system_torque;
	extern torque_statistic_t grain_torque;
	extern torque_statistic_t material_torque;

   extern specific_heat_statistic_t system_specific_heat;
	extern specific_heat_statistic_t grain_specific_heat;
   extern specific_heat_statistic_t material_specific_heat;

   extern susceptibility_statistic_t system_susceptibility;
	extern susceptibility_statistic_t grain_susceptibility;
   extern susceptibility_statistic_t material_susceptibility;

   extern standard_deviation_statistic_t material_standard_deviation;

   extern spin_length_statistic_t system_spin_length;
   extern spin_length_statistic_t material_spin_length;
   extern spin_length_statistic_t height_spin_length;

   extern binder_cumulant_statistic_t system_binder_cumulant;
   extern binder_cumulant_statistic_t material_binder_cumulant;

}

#endif /*STATS_H_*/
