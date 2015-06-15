#ifndef CUDA_INTERNAL_H_
#define CUDA_INTERNAL_H_
//-----------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------

//---------------------------------------------------------------------
// Defines shared internal data structures and functions for the
// cuda implementation. These functions should
// not be accessed outside of the local temperature pulse code.
//---------------------------------------------------------------------

#include <thrust/copy.h>
#include <thrust/device_vector.h>

/*
 * requesting data strcutures from the main program
 */

#include "atoms.hpp"
#include "cells.hpp"
#include "material.hpp"

namespace cuda{
   namespace internal{

      //-----------------------------------------------------------------------------
      // Shared variables used for the cuda implementation
      //-----------------------------------------------------------------------------

      /*
       * Atom information
       */
      thrust::device_vector<double> x_spin_array;
      thrust::device_vector<double> y_spin_array;
      thrust::device_vector<double> z_spin_array;

      thrust::device_vector<double> x_coord_array;
      thrust::device_vector<double> y_coord_array;
      thrust::device_vector<double> z_coord_array;

      thrust::device_vector<size_t> type_array;

      thrust::device_vector<size_t> cell_array;

      bool __initialize_atoms ();

      /*
       * Field information
       */

      thrust::device_vector<double> x_total_spin_field_array;
      thrust::device_vector<double> y_total_spin_field_array;
      thrust::device_vector<double> z_total_spin_field_array;

      thrust::device_vector<double> x_total_external_field_array;
      thrust::device_vector<double> y_total_external_field_array;
      thrust::device_vector<double> z_total_external_field_array;

      /*
       * Required by the total external field calculator
       * and the dipolar field updater
       */
      thrust::device_vector<double> x_dipolar_field_array;
      thrust::device_vector<double> y_dipolar_field_array;
      thrust::device_vector<double> z_dipolar_field_array;

      bool __initialize_fields ();

      /*
       * Cell information
       */

      thrust::device_vector<double> cell_x_coord_array;
      thrust::device_vector<double> cell_y_coord_array;
      thrust::device_vector<double> cell_z_coord_array;

      thrust::device_vector<double> cell_x_mag_array;
      thrust::device_vector<double> cell_y_mag_array;
      thrust::device_vector<double> cell_z_mag_array;

      thrust::device_vector<double> cell_volume_array;

      thrust::device_vector<size_t> cell_num_atoms;

      bool __initialize_cells ();

      /*
       * Material information
       */

      thrust::device_vector<mp::materials_t> materials;

      bool __initialize_materials ();

      /*
       * Topology information
       */

      thrust::device_vector<size_t> limits;
      thrust::device_vector<size_t> neighbours;

      bool __initialize_topology ();

      /*
       * Functors
       */

      struct plusone_functor
      {
         __host__ __device__
            size_t operator() (const float& item) const
            {
               return item + 1UL;
            }
      };

      //-----------------------------------------------------------------------------
      // Shared functions and kernels used for the cuda implementation
      //-----------------------------------------------------------------------------

      __global__ void llg_heun_first_kernel (
            double * x_spin, double * y_spin, double * z_spin,
            double * x_sp_field, double * y_sp_field, double * z_sp_field,
            double * x_ext_field, double * y_ext_field, double * z_ext_field,
            double dt
            );
     __global__ void llg_heun_scheme(
	    double * x_spin, double * y_spin, double * z_spin,   
	    double * x_sp_field, double * y_sp_field, double * z_sp_field,
	    double * x_ext_field, double * y_ext_field, double * z_ext_field,
	    double * x_new_spin, double * y_new_spin, double z_new_spin	
            );
   } // end of iternal namespace
} // end of cuda namespace

#endif //CUDA_INTERNAL_H_
