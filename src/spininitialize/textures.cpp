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

// C++ standard library headers
#include <algorithm>
#include <cmath>

// Vampire headers
#include "spininitialize.hpp"
#include "random.hpp"

// spininitialize module headers
#include "internal.hpp"

namespace spininitialize{

   namespace internal{

      namespace {
         const double pi = 3.14159265358979323846;
      }

      //-----------------------------------------------------------------------------
      // uniform_vector: spins point along the user-defined direction mat.initial_spin
      //-----------------------------------------------------------------------------
      void uniform_vector_spin(const mp_t& mat, double& sx, double& sy, double& sz){

         sx = mat.initial_spin[0];
         sy = mat.initial_spin[1];
         sz = mat.initial_spin[2];

      }

      //-----------------------------------------------------------------------------
      // random: spins point in a random direction, uniformly distributed on the unit
      // sphere (via normalised gaussian components)
      //-----------------------------------------------------------------------------
      void random_spin(MTRand& prng, double& sx, double& sy, double& sz){

         sx = mtrandom::gaussianc(prng);
         sy = mtrandom::gaussianc(prng);
         sz = mtrandom::gaussianc(prng);

      }

      //-----------------------------------------------------------------------------
      // domain_wall: 1D tanh-profile Neel-type domain wall.
      //
      // Along the chosen axis the spin rotates from +<axis> (far below the wall
      // centre) through the rotation axis itself (at the wall centre) to
      // -<axis> (far above the wall centre). The rotation takes place in the
      // plane spanned by <axis> and the next axis in cyclic order (x->y->z->x),
      // i.e. for domain-wall-x the rotation plane is x-z, for domain-wall-y it
      // is y-x and for domain-wall-z it is z-y.
      //
      // mat.centre[0] is the position of the wall centre (fraction of system
      // size along axis) and mat.width controls the sharpness of the wall
      // (fraction of system size); smaller values give a sharper wall.
      //-----------------------------------------------------------------------------
      void domain_wall_spin(const mp_t& mat, const double fx, const double fy, const double fz, double& sx, double& sy, double& sz){

         const double pos[3] = {fx, fy, fz};

         const double t = (pos[mat.axis] - mat.centre[0]) / mat.width;

         // theta varies smoothly from 0 (t -> -inf) to pi (t -> +inf)
         const double theta = 2.0 * atan(exp(t));

         double m[3] = {0.0, 0.0, 0.0};
         m[mat.axis]       = cos(theta); // component along the wall axis
         m[(mat.axis+2)%3] = sin(theta); // in-plane rotation component

         sx = m[0];
         sy = m[1];
         sz = m[2];

      }

      //-----------------------------------------------------------------------------
      // skyrmion: radially symmetric "cone" skyrmion in the x-y plane.
      //
      // mat.centre = (cx, cy) is the centre of the skyrmion (fraction of system
      // size), mat.width is its radius (fraction of system size). mat.polarity
      // (+/-1) sets the out-of-plane direction of the core, and mat.chirality
      // (+/-1) selects between a Neel-type (radial, +1) and Bloch-type
      // (azimuthal, -1) in-plane rotation sense.
      //-----------------------------------------------------------------------------
      void skyrmion_spin(const mp_t& mat, const double fx, const double fy, const double fz, double& sx, double& sy, double& sz){

         const double dx = fx - mat.centre[0];
         const double dy = fy - mat.centre[1];

         const double r = sqrt(dx*dx + dy*dy);
         const double phi = atan2(dy, dx);

         const double R = mat.width;

         // theta(r) varies linearly from pi at the core (r=0) to 0 at the
         // skyrmion edge (r=R) and beyond
         const double theta = pi * (1.0 - std::min(r / R, 1.0));

         // Neel (chirality +1) or Bloch (chirality -1) in-plane rotation
         const double gamma = (mat.chirality >= 0) ? 0.0 : (0.5 * pi);

         sx = sin(theta) * cos(phi + gamma);
         sy = sin(theta) * sin(phi + gamma);
         sz = -double(mat.polarity) * cos(theta);

      }

      //-----------------------------------------------------------------------------
      // spin_spiral: helical spin spiral propagating along the chosen axis.
      //
      // The spin rotates in the plane spanned by the two axes perpendicular to
      // the propagation axis (cyclically: x propagation -> y-z rotation plane,
      // y -> z-x, z -> x-y), completing one full turn over mat.width (the
      // wavelength, as a fraction of the system size along the axis).
      //-----------------------------------------------------------------------------
      void spin_spiral_spin(const mp_t& mat, const double fx, const double fy, const double fz, double& sx, double& sy, double& sz){

         const double pos[3] = {fx, fy, fz};

         const double k = 2.0 * pi * pos[mat.axis] / mat.width;

         double m[3] = {0.0, 0.0, 0.0};
         m[mat.axis]       = 0.0;
         m[(mat.axis+1)%3] = cos(k);
         m[(mat.axis+2)%3] = sin(k);

         sx = m[0];
         sy = m[1];
         sz = m[2];

      }

      //-----------------------------------------------------------------------------
      // vortex: in-plane curling vortex in the x-y plane, with an optional
      // out-of-plane core.
      //
      // mat.centre = (cx, cy) is the centre of the vortex (fraction of system
      // size). mat.width is the core radius (fraction of system size); the
      // out-of-plane component of the core decays as a gaussian with this
      // radius and vanishes for mat.width = 0 (a pure in-plane vortex).
      // mat.polarity (+/-1) sets the out-of-plane core direction and
      // mat.chirality (+/-1) sets the sense (clockwise/anticlockwise) of the
      // in-plane curling.
      //-----------------------------------------------------------------------------
      void vortex_spin(const mp_t& mat, const double fx, const double fy, const double fz, double& sx, double& sy, double& sz){

         const double dx = fx - mat.centre[0];
         const double dy = fy - mat.centre[1];

         const double r = sqrt(dx*dx + dy*dy);
         const double phi = atan2(dy, dx);

         double mz = 0.0;
         if(mat.width > 0.0){
            const double u = r / mat.width;
            mz = double(mat.polarity) * exp(-u*u);
         }

         const double in_plane = sqrt(std::max(0.0, 1.0 - mz*mz));

         sx = -sin(phi) * double(mat.chirality) * in_plane;
         sy =  cos(phi) * double(mat.chirality) * in_plane;
         sz = mz;

      }

      //-----------------------------------------------------------------------------
      // Master dispatch function returning the (unnormalised) initial spin
      // direction for material mat at fractional coordinates (fx,fy,fz).
      //-----------------------------------------------------------------------------
      void get_spin_direction(const int mat, const double fx, const double fy, const double fz,
                               double& sx, double& sy, double& sz, MTRand& prng){

         // default direction for materials with no spininitialize entry
         if(mat < 0 || mat >= int(mp.size())){
            sx = 0.0;
            sy = 0.0;
            sz = 1.0;
            return;
         }

         const mp_t& material = mp[mat];

         switch(material.texture){

            case random:
               random_spin(prng, sx, sy, sz);
               break;

            case domain_wall:
               domain_wall_spin(material, fx, fy, fz, sx, sy, sz);
               break;

            case skyrmion:
               skyrmion_spin(material, fx, fy, fz, sx, sy, sz);
               break;

            case spin_spiral:
               spin_spiral_spin(material, fx, fy, fz, sx, sy, sz);
               break;

            case vortex:
               vortex_spin(material, fx, fy, fz, sx, sy, sz);
               break;

            case vector_field:
               vector_field_spin(material, fx, fy, fz, sx, sy, sz);
               break;

            case uniform_vector:
            default:
               uniform_vector_spin(material, sx, sy, sz);
               break;

         }

         // normalise the resulting direction vector, falling back to +z if it is
         // degenerate (e.g. exactly at a singular point of a texture)
         if(!normalise(sx, sy, sz)){
            sx = 0.0;
            sy = 0.0;
            sz = 1.0;
         }

      }

   } // end of internal namespace

} // end of spininitialize namespace
