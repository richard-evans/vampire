//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef SPININITIALIZE_INTERNAL_H_
#define SPININITIALIZE_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the spininitialize module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <string>
#include <vector>

// Vampire headers
#include "spininitialize.hpp"
#include "random.hpp"

// spininitialize module headers
#include "internal.hpp"

namespace spininitialize{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      //-----------------------------------------------------------------------------
      // Enumeration of the different supported initial spin textures.
      //
      // Each material has exactly one texture, selected via
      // material[#]:initial-spin-direction in the material file (see
      // interface.cpp). The texture value determines which of the
      // "individual texture functions" declared below (and implemented in
      // textures.cpp / vector_field.cpp) is used to compute the initial
      // spin direction of each atom of that material.
      //-----------------------------------------------------------------------------
      enum spin_texture_t{
         uniform_vector = 0, // single user-defined unit vector (default (0,0,1))
         random         = 1, // spins initialised to random directions (infinite temperature)
         domain_wall    = 2, // 1D tanh-profile domain wall along x, y or z
         skyrmion       = 3, // radially symmetric skyrmion texture in the x-y plane
         spin_spiral    = 4, // helical spin spiral propagating along x, y or z
         vortex         = 5, // in-plane curling vortex texture in the x-y plane
         vector_field   = 6  // texture interpolated from a user-supplied vector field file
      };

      //-----------------------------------------------------------------------------
      // A single point of a user-defined vector field
      // (x,y,z) are fractional coordinates of the system size [0:1]
      // (mx,my,mz) is the (not necessarily normalised) spin direction at that point
      //-----------------------------------------------------------------------------
      class field_point_t{
         public:
            double x,y,z;
            double mx,my,mz;
            field_point_t(): x(0.0), y(0.0), z(0.0), mx(0.0), my(0.0), mz(1.0) {}
      };

      //-----------------------------------------------------------------------------
      // internal materials class for storing material spin initialisation parameters
      //-----------------------------------------------------------------------------
      class mp_t{

          private:

          public:

             //------------------------------
             // material parameter variables
             //------------------------------

             // type of spin texture for this material
             spin_texture_t texture;

             // uniform_vector: direction of the spins
             double initial_spin[3];

             // domain_wall: axis along which the wall varies (0=x, 1=y, 2=z)
             // spin_spiral: axis along which the spiral propagates (0=x, 1=y, 2=z)
             int axis;

             // domain_wall: centre[0] = position of the wall centre (fraction of system size along axis)
             // skyrmion, vortex: centre[0],centre[1] = centre of the texture in the x-y plane (fraction of system size)
             double centre[2];

             // domain_wall: width of the wall (fraction of system size along axis)
             // skyrmion, vortex: radius of the texture (fraction of system size)
             // spin_spiral: wavelength of the spiral (fraction of system size along axis)
             double width;

             // skyrmion, vortex: sense of in-plane rotation (+1 or -1)
             int chirality;

             // skyrmion, vortex: out-of-plane direction of the core (+1 or -1)
             int polarity;

             // vector_field: index into internal::vector_field_data / internal::vector_field_filenames
             int vector_field_id;

             // constructor
             mp_t (const unsigned int max_materials = 100):
                texture(uniform_vector),
                axis(0),
                width(0.1),
                chirality(1),
                polarity(1),
                vector_field_id(-1)
             {
                initial_spin[0] = 0.0;
                initial_spin[1] = 0.0;
                initial_spin[2] = 1.0;
                centre[0] = 0.5;
                centre[1] = 0.5;
             }; // end of constructor

       }; // end of internal::mp class

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------

      extern bool enabled; // bool to enable module

      extern std::vector<internal::mp_t> mp; // array of material spin initialisation properties

      // cache of loaded vector field files, indexed by internal::mp[mat].vector_field_id
      extern std::vector<std::string> vector_field_filenames;
      extern std::vector< std::vector<field_point_t> > vector_field_data;

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

      // split a comma-separated list of values into trimmed string tokens
      std::vector<std::string> split_csv(const std::string& value);

      // load (or retrieve from cache) a vector field file, returning its id
      int get_vector_field_id(const std::string& filename, const int line);

      // normalise a 3-vector in place; returns false if the vector has zero length
      bool normalise(double& sx, double& sy, double& sz);

      // dispatch function returning the initial spin direction (unnormalised) for a given material and fractional position
      void get_spin_direction(const int mat, const double fx, const double fy, const double fz,
                               double& sx, double& sy, double& sz, MTRand& prng);

      // individual texture functions, all working in fractional coordinates [0:1]
      void uniform_vector_spin(const mp_t& mat, double& sx, double& sy, double& sz);
      void random_spin(MTRand& prng, double& sx, double& sy, double& sz);
      void domain_wall_spin(const mp_t& mat, const double fx, const double fy, const double fz, double& sx, double& sy, double& sz);
      void skyrmion_spin(const mp_t& mat, const double fx, const double fy, const double fz, double& sx, double& sy, double& sz);
      void spin_spiral_spin(const mp_t& mat, const double fx, const double fy, const double fz, double& sx, double& sy, double& sz);
      void vortex_spin(const mp_t& mat, const double fx, const double fy, const double fz, double& sx, double& sy, double& sz);
      void vector_field_spin(const mp_t& mat, const double fx, const double fy, const double fz, double& sx, double& sy, double& sz);

   } // end of internal namespace

} // end of spininitialize namespace

#endif //SPININITIALIZE_INTERNAL_H_
