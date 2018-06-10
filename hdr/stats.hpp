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
	extern int num_atoms;				/// Number of atoms for statistic purposes
	extern double inv_num_atoms;	///1.0/num_atoms
	extern double max_moment;		/// Total Maximum moment
	extern double data_counter;		/// number of data points for averaging

	// Member Functions
	extern int mag_m();
	extern void mag_m_reset();
	extern double max_torque();

	extern bool calculate_torque;
	extern double total_system_torque[3];
	extern double total_mean_system_torque[3];

	extern std::vector <double> sublattice_mean_torque_x_array;
	extern std::vector <double> sublattice_mean_torque_y_array;
	extern std::vector <double> sublattice_mean_torque_z_array;

	extern double torque_data_counter;

   extern bool calculate_energy;

   /// Statistics energy types
   enum energy_t { total = 0, exchange = 1, anisotropy = 2, applied_field = 3, magnetostatic = 4};

   /// Statistics types
   enum stat_t { atotal=0, mean=1};

   /// Statistics output functions
   extern void output_energy(std::ostream&, enum energy_t, enum stat_t);

   //-------------------------------------------------
   // New statistics module functions and variables
   //-------------------------------------------------

   // Control functions
   void initialize(const int num_atoms, const int num_materials,
                   const std::vector<double>& magnetic_moment_array,
                   const std::vector<int>& material_type_array,
                   const std::vector<int>& height_category_array,
                   const std::vector<bool>& non_magnetic_materials_array);

   void update(const std::vector<double>& sx, const std::vector<double>& sy, const std::vector<double>& sz,
               const std::vector<double>& mm, const std::vector<int>& mat, const double temperature);
   void reset();

   // Statistics control flags (to be moved internally when long-awaited refactoring of vio is done)
   extern bool calculate_system_energy;
   extern bool calculate_system_magnetization;
   extern bool calculate_material_magnetization;
   extern bool calculate_height_magnetization;
   extern bool calculate_material_height_magnetization;
   extern bool calculate_system_susceptibility;
   extern bool calculate_material_susceptibility;

   // forward declaration of friend classes
   class susceptibility_statistic_t;

   //----------------------------------
   // Energy class definition
   //----------------------------------
   class energy_statistic_t{

      //friend class specific_heat_statistic_t;

   public:
      energy_statistic_t ();
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

      std::string output_energy(enum energy_t energy_type);
      std::string output_mean_energy(enum energy_t energy_type);

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

   };

   //----------------------------------
   // Magnetization Class definition
   //----------------------------------
   class magnetization_statistic_t{

      friend class susceptibility_statistic_t;

      public:
         magnetization_statistic_t ();
         bool is_initialized();
         void set_mask(const int mask_size, std::vector<int> inmask, const std::vector<double>& mm);
         void get_mask(std::vector<int>& out_mask, std::vector<double>& out_saturation);
         void calculate_magnetization(const std::vector<double>& sx, const std::vector<double>& sy, const std::vector<double>& sz, const std::vector<double>& mm);
         void set_magnetization(std::vector<double>& magnetization, std::vector<double>& mean_magnetization, long counter);
         void reset_magnetization_averages();
         const std::vector<double>& get_magnetization();
         std::string output_magnetization();
         std::string output_normalized_magnetization();
         std::string output_normalized_magnetization_length();
         std::string output_normalized_mean_magnetization();
         std::string output_normalized_mean_magnetization_length();
         std::string output_normalized_magnetization_dot_product(const std::vector<double>& vec);
         std::string output_mean_magnetization_length();
         std::string output_mean_magnetization();

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

   };

   //----------------------------------
   // Susceptibility Class definition
   //----------------------------------
   class susceptibility_statistic_t{

      public:
         susceptibility_statistic_t ();
         void initialize(magnetization_statistic_t& mag_stat);
         void calculate(const std::vector<double>& magnetization);
         void reset_averages();
         std::string output_mean_susceptibility(const double temperature);
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

   };

   // Statistics classes
   extern energy_statistic_t system_energy;
   extern energy_statistic_t material_energy;

   extern magnetization_statistic_t system_magnetization;
   extern magnetization_statistic_t material_magnetization;
   extern magnetization_statistic_t height_magnetization;
   extern magnetization_statistic_t material_height_magnetization;

   extern susceptibility_statistic_t system_susceptibility;
   extern susceptibility_statistic_t material_susceptibility;

}

#endif /*STATS_H_*/
