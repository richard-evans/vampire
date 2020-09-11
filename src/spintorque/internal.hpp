#ifndef SPINTORQUE_INTERNAL_H_
#define SPINTORQUE_INTERNAL_H_
//-----------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans and P Chureemart 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

//---------------------------------------------------------------------
// Defines shared internal data structures and functions for the
// spin-torque implementation. These functions should not be accessed
// outside of the spin-torque code.
//---------------------------------------------------------------------
namespace st{
   namespace internal{

      //-----------------------------------------------------------------------------
      // Shared variables used for the spin torque calculation
      //-----------------------------------------------------------------------------
      extern bool enabled; // enable spin torque calculation
      extern bool TMRenable;

      extern double micro_cell_size; /// lateral size of spin torque microcells
      extern double micro_cell_thickness; /// thickness of spin torque microcells (atomistic)

      extern int num_local_atoms; /// number of local atoms (ignores halo atoms in parallel simulation)
      extern int current_direction; /// direction for current x->0, y->1, z->2
      //   std::vector< std::vector< micro_cell_t > > stack;
      extern std::vector<int> atom_st_index; // mc which atom belongs to
      extern std::vector<double> x_field_array; // arrays to store atomic spin torque field
      extern std::vector<double> y_field_array;
      extern std::vector<double> z_field_array;

      extern int num_stacks;  // total number of stacks
      extern int num_x_stacks; // number of stacks in x
      extern int num_y_stacks; // number of stack in y
      extern int num_microcells_per_stack; // number of microcells per stack

      extern int config_file_counter; // spin torque config file counter

      extern int free_layer;       /// index of free layer in magnetic tunnel junction
      extern int reference_layer;  /// index of reference layer in magnetic tunnel junction

      extern double je; // current (C/s)
      extern double initial_beta;
      extern double rel_angle;
      extern int ST_output_rate;

      extern std::vector<double> initial_m;


      extern std::vector<int> stack_index; // start of stack in microcell arrays

      extern std::vector<double> beta_cond; /// spin polarisation (conductivity)
      extern std::vector<double> beta_diff; /// spin polarisation (diffusion)
      extern std::vector<double> sa_infinity; /// intrinsic spin accumulation
      extern std::vector<double> lambda_sdl; /// spin diffusion length
      extern std::vector<double> diffusion; /// spin diffusion length
      extern std::vector<double> sd_exchange; /// spin diffusion length
      extern std::vector<double> a; /// spin diffusion length
      extern std::vector<double> b; /// spin diffusion length
      extern std::vector<double> coeff_ast;
      extern std::vector<double> coeff_nast;

      extern std::vector<double> cell_natom;


      // three-vector arrays
      extern std::vector<double> pos; /// stack position
      extern std::vector<double> m; // magnetisation
      extern std::vector<double> j; // spin current
      extern std::vector<double> sa; // spin accumulation
      extern std::vector<double> spin_torque; // spin torque
      extern std::vector<double> ast; // adiabatic spin torque
      extern std::vector<double> nast; // non-adiabatic spin torque
      extern std::vector<double> total_ST; // non-adiabatic spin torque
      extern std::vector<double> magx_mat; // magnetisation of material
      extern std::vector<double> magy_mat;
      extern std::vector<double> magz_mat;
      
      

      // material parameters for spin torque calculation
      struct mp_t{
         double beta_cond;    /// spin polarisation (conductivity)
         double beta_diff;    /// spin polarisation (diffusion)
         double sa_infinity;  /// intrinsic spin accumulation
         double lambda_sdl;   /// spin diffusion length
         double diffusion;    /// diffusion constant
         double sd_exchange;  /// sd_exchange constant
      };

      // three vector type definition
      class three_vector_t{
      public:
         double x;
         double y;
         double z;

         // constructor
         three_vector_t(double ix, double iy, double iz){
            x = ix;
            y = iy;
            z = iz;
         }

      };

      // matrix type definition
      struct matrix_t{
         double xx;
         double xy;
         double xz;
         double yx;
         double yy;
         double yz;
         double zx;
         double zy;
         double zz;
      };

      // array of material properties
      extern std::vector<st::internal::mp_t> mp;

      // default material properties
      extern st::internal::mp_t default_properties;

      //-----------------------------------------------------------------------------
      // Shared functions used for the spin torque calculation
      //-----------------------------------------------------------------------------
      void output_microcell_data();
      void output_base_microcell_data();
      void calculate_spin_accumulation();
      void update_cell_magnetisation(const std::vector<double>& x_spin_array,
                                     const std::vector<double>& y_spin_array,
                                     const std::vector<double>& z_spin_array,
                                     const std::vector<int>& atom_type_array,
                                     const std::vector<double>& mu_s_array);

      void set_inverse_transformation_matrix(const st::internal::three_vector_t& reference_vector, st::internal::matrix_t& itm);
      st::internal::three_vector_t transform_vector(const st::internal::three_vector_t& rv, const st::internal::matrix_t& tm);
      st::internal::three_vector_t gaussian_elimination(st::internal::matrix_t& M, st::internal::three_vector_t& V);


   } // end of iternal namespace
} // end of st namespace

#endif //SPINTORQUE_INTERNAL_H_
