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

         // every atom of this material gets the same fixed direction, as
         // entered by the user (e.g. "0,0,1" or "[110]"); this has already
         // been normalised to a unit vector by check_for_valid_unit_vector
         // when the input file was read
         sx = mat.initial_spin[0];
         sy = mat.initial_spin[1];
         sz = mat.initial_spin[2];

      }

      //-----------------------------------------------------------------------------
      // random: spins point in a random direction, uniformly distributed on the unit
      // sphere (via normalised gaussian components)
      //-----------------------------------------------------------------------------
      void random_spin(MTRand& prng, double& sx, double& sy, double& sz){

         // drawing each Cartesian component independently from a gaussian
         // distribution and then normalising the resulting vector (done by
         // the caller, get_spin_direction) gives a direction that is
         // uniformly distributed over the surface of the unit sphere. This
         // is NOT the same as drawing sx,sy,sz uniformly from [-1,1], which
         // would bias spins towards the corners of the cube.
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

         // pick out the fractional coordinate along the wall's propagation
         // axis (0=x, 1=y, 2=z), as set by domain-wall-x/y/z
         const double pos[3] = {fx, fy, fz};

         // t is a dimensionless "distance from the wall centre" measured in
         // units of the wall width: t = 0 at the centre, |t| >> 1 far away.
         // mat.centre[0] and mat.width are both fractions of the system size.
         const double t = (pos[mat.axis] - mat.centre[0]) / mat.width;

         // standard tanh-profile domain wall: the polar angle theta of the
         // spin (measured from the +axis direction) follows
         //   theta(t) = 2*atan(exp(t))
         // which smoothly varies from 0 (spin along +axis, t -> -inf)
         // through pi/2 (spin perpendicular to axis, t = 0, the wall centre)
         // to pi (spin along -axis, t -> +inf).
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

         // (dx,dy) is the in-plane displacement of this atom from the
         // skyrmion centre (cx,cy), in fractional coordinates
         const double dx = fx - mat.centre[0];
         const double dy = fy - mat.centre[1];

         // r is the radial distance from the centre, phi is the azimuthal
         // (in-plane) angle around the centre
         const double r = sqrt(dx*dx + dy*dy);
         const double phi = atan2(dy, dx);

         const double R = mat.width; // skyrmion radius (fraction of system size)

         // theta(r) is the polar angle of the spin away from the z axis.
         // It varies linearly from pi at the core (r=0, spin antiparallel to
         // the core direction set below) to 0 at the skyrmion edge (r=R)
         // and beyond, where spins simply point along +/-z (the background).
         // std::min(r/R, 1.0) clamps theta to 0 outside the skyrmion radius.
         const double theta = pi * (1.0 - std::min(r / R, 1.0));

         // gamma is an extra in-plane rotation applied to the azimuthal
         // angle phi, which sets the "type" of skyrmion:
         //   chirality = +1 -> gamma = 0      : Neel-type (radial in-plane spins)
         //   chirality = -1 -> gamma = pi/2   : Bloch-type (azimuthal in-plane spins)
         const double gamma = (mat.chirality >= 0) ? 0.0 : (0.5 * pi);

         // convert (theta, phi+gamma) spherical angles to a Cartesian unit
         // vector. At the core (theta=pi) sz = +polarity, i.e. the core
         // points along the chosen polarity direction; far from the core
         // (theta=0) sz = -polarity, the opposite (background) direction.
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

         // fractional coordinate along the propagation axis (0=x, 1=y, 2=z)
         const double pos[3] = {fx, fy, fz};

         // k is the phase of the spiral at this position: it increases by
         // 2*pi (one full rotation of the spin) every time pos[axis]
         // advances by mat.width (the wavelength, as a fraction of the
         // system size along the propagation axis)
         const double k = 2.0 * pi * pos[mat.axis] / mat.width;

         // the spin has no component along the propagation axis; instead it
         // rotates within the plane perpendicular to it, tracing out a helix
         // as a function of position
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

         // (dx,dy) is the in-plane displacement of this atom from the
         // vortex centre (cx,cy), in fractional coordinates
         const double dx = fx - mat.centre[0];
         const double dy = fy - mat.centre[1];

         // r is the radial distance from the centre, phi is the azimuthal
         // (in-plane) angle around the centre
         const double r = sqrt(dx*dx + dy*dy);
         const double phi = atan2(dy, dx);

         // out-of-plane component mz: zero everywhere for a pure in-plane
         // vortex (mat.width = 0). For mat.width > 0, mz follows a gaussian
         // that peaks at the vortex centre (mz = +/-1 at r=0, set by
         // mat.polarity) and decays to zero over a radius ~ mat.width,
         // modelling the localised out-of-plane core seen in real vortices.
         double mz = 0.0;
         if(mat.width > 0.0){
            const double u = r / mat.width;
            mz = double(mat.polarity) * exp(-u*u);
         }

         // the remaining in-plane magnitude, so that (in_plane)^2 + mz^2 = 1
         const double in_plane = sqrt(std::max(0.0, 1.0 - mz*mz));

         // the in-plane component curls azimuthally around the centre
         // (tangential to circles of constant r), i.e. perpendicular to the
         // radial direction (cos(phi), sin(phi)). mat.chirality = +1 gives
         // anticlockwise curling, -1 gives clockwise.
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
         // (e.g. if the mp array was never resized for this material because
         // it has no initial-spin-direction keyword in the material file)
         if(mat < 0 || mat >= int(mp.size())){
            sx = 0.0;
            sy = 0.0;
            sz = 1.0;
            return;
         }

         const mp_t& material = mp[mat];

         // dispatch to the appropriate texture function based on the
         // texture type selected for this material when the input file was
         // parsed (see interface.cpp)
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
