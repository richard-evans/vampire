//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <chrono>
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <vector>

// Vampire headers
#include "atoms.hpp"
#include "cells.hpp"
#include "errors.hpp"
#include "gpu.hpp"
#include "material.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vopencl.hpp"

// vopencl module headers
#include "cl/material_type.h"
#include "data.hpp"
#include "init_kernels.hpp"
#include "internal.hpp"
#include "opencl_include.hpp"
#include "opencl_utils.hpp"
#include "statistics.hpp"
#include "typedefs.hpp"

#ifdef OPENCL
namespace vcl = ::vopencl::internal;

namespace vopencl
{
   namespace internal
   {
      cl::NDRange global;
      cl::NDRange global_other;

      namespace time
      {
         time_t sim_start;

#ifdef OPENCL_TIME_KERNELS
         double spin_fields = 0.0;
         double mat_mul = 0.0;
         double rng = 0.0;
         double external_fields = 0.0;
         double predictor_step = 0.0;
         double corrector_step = 0.0;
#endif // OPENCL_TIME_KERNELS
      }
   }
}
#endif // OPENCL

namespace vopencl
{
   //----------------------------------------------------------------------------
   // Function to initialize vopencl module
   //----------------------------------------------------------------------------
   bool initialize(bool cpu_stats)
   {
      bool success = false;

#ifdef OPENCL

      // initialisation start time
      auto start = std::chrono::high_resolution_clock::now();

      std::string message("OpenCL has been enabled in "
#ifdef OPENCL_DP
                          "double precision mode."
#else
                          "single precision mode."
#endif // OPENCL_DP
 
#ifdef OPENCL_USE_VECTOR_TYPE
                          " OpenCL vector types will be used for storage."
#endif // OPENCL_USE_VECTOR_TYPE

#ifdef OPENCL_USE_NATIVE_FUNCTIONS
                          " Native functions will be used."
#endif // OPENCL_USE_NATIVE_FUNCTIONS

#ifdef OPENCL_DEBUG
                          " Debugging routines have been enabled."
#endif // OPENCL_DEBUG

#ifdef OPENCL_TIME_KERNELS
                          " OpenCL kernels will be timed."
#endif // OPENCL_TIME_KERNELS
         );

      std::cout << message << std::endl;
      zlog << zTs() << message << std::endl;

      //vcl::stats::use_cpu = cpu_stats;
      vcl::stats::use_cpu = true;

      // find OpenCL platforms
      std::vector<cl::Platform> platforms;
      cl::Platform::get(&platforms);
      unsigned nplatforms = platforms.size();

      if (nplatforms == 0)
      {
         message = "Error: OpenCL is enabled but no platforms are available.";
         terminaltextcolor(RED);
         std::cout << message << std::endl;
         terminaltextcolor(WHITE);
         zlog << zTs() << message << std::endl;
         ::err::vexit();
      }

      // find available devices for each platform
      std::vector<std::vector<cl::Device>> devices(nplatforms);
      unsigned ndevices = 0;
      for (unsigned i=0; i<nplatforms; ++i)
      {
         std::vector<cl::Device> tmp_devices;
         platforms[i].getDevices(CL_DEVICE_TYPE_ALL, &tmp_devices);
         devices[i] = tmp_devices;
         ndevices += tmp_devices.size();

#ifdef OPENCL_LOG
         vcl::OCLLOG << "Found platform " << platforms[i].getInfo<CL_PLATFORM_NAME>() << std::endl;
         for (unsigned j=0; j<tmp_devices.size(); ++j)
         {
            vcl::OCLLOG << "Found device " << tmp_devices[j].getInfo<CL_DEVICE_NAME>() << std::endl;
            vcl::OCLLOG << "with version " << tmp_devices[j].getInfo<CL_DEVICE_VERSION>() << std::endl;
         }
#endif // OPENCL_LOG
      }

      if (ndevices == 0)
      {
         message = "Error: OpenCL is enabled but no suitable devices can be found.";
         terminaltextcolor(RED);
         std::cout << message << std::endl;
         terminaltextcolor(WHITE);
         zlog << zTs() << message << std::endl;
         ::err::vexit();
      }

      //// Set up platforms and devices ////
      if (::gpu::platform < 0) // not specified in input file
      {
         // use first platform found
         ::gpu::platform = 0;
      }
      else if (::gpu::platform >= nplatforms) // specified incorrectly
      {
         terminaltextcolor(YELLOW);
         std::cerr << "Warning: Platform specified does not exist (" << ::gpu::platform << "). "
                   << "Falling back to platform 0." << std::endl;
         terminaltextcolor(WHITE);
         ::gpu::platform = 0;
      }

      if (::gpu::device < 0) // not specified in input file
      {
         // silently use default
         ::gpu::device = 0;
      }
      else if (::gpu::device >= devices[::gpu::platform].size()) // specified incorrectly
      {
         terminaltextcolor(YELLOW);
         std::cerr << "Warning: Device specified (" << ::gpu::device
                   << ") does not exist on platform " << ::gpu::platform << ". "
                   << "Falling back to device 0." << std::endl;
         terminaltextcolor(WHITE);
         ::gpu::device = 0;
      }

      if (::gpu::num_threads > 0)
      {
         vcl::global = cl::NDRange(::gpu::num_threads);
      }
      else
      {
         int default_nthreads = 8192;
         terminaltextcolor(YELLOW);
         std::cerr << "Warning: defaulting to " << default_nthreads << " threads on main device." << std::endl;
         std::cerr << "Use gpu:num-threads=n in the input file to specify number." << std::endl;
         terminaltextcolor(WHITE);
         vcl::global = cl::NDRange(default_nthreads);
      }

      // platform and device to use for the majority of kernels
      cl::Platform default_platform = platforms[::gpu::platform];
      vcl::default_device = devices[::gpu::platform][::gpu::device];

#ifndef ENABLE_MULTIPLE_DEVICES
      // if not using multiple devices for kernels
      // i.e. just a host CPU and a device, which may be the same CPU

      vcl::context = cl::Context({vcl::default_device});

      vcl::queue = cl::CommandQueue(vcl::context, vcl::default_device,
                                    CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE);
#else
      // using multiple OpenCL devices

      // which platform to use for extra device
      if (::gpu::platform_other < 0) // not specified in input file
      {
         // use the platform already being used
         ::gpu::platform_other = ::gpu::platform;
      }
      else if (::gpu::platform_other >= nplatforms) // specified incorrectly
      {
         terminaltextcolor(YELLOW);
         std::cerr << "Warning: Platform specified does not exist (" << ::gpu::platform << "). "
                   << "Falling back to platform " << ::gpu::platform << '.' << std::endl;
         terminaltextcolor(WHITE);
         ::gpu::platform_other = ::gpu::platform;
      }

      // which device to use
      if (::gpu::device_other < 0) // not specified in input file
      {
         // silently use only one device
         // the same as if ENABLE_MULTIPLE_DEVICES wasn't defined
         ::gpu::device_other = ::gpu::device;
      }
      else if (::gpu::device_other >= devices[::gpu::platform_other].size()) // specified incorrectly
      {
         terminaltextcolor(YELLOW);
         std::cerr << "Warning: Device specified (" << ::gpu::device_other
                   << ") does not exist on platform " << ::gpu::platform_other << ". "
                   << "Falling back to device " << ::gpu::device << '.' << std::endl;
         terminaltextcolor(WHITE);
         ::gpu::device_other = ::gpu::device;
      }

      if (::gpu::num_threads_other > 0)
      {
         vcl::global_other = cl::NDRange(::gpu::num_threads_other);
      }
      else
      {
         vcl::global_other = vcl::global;
      }

      cl::Platform extra_platform = platforms[::gpu::platform_other];
      vcl::extra_device = devices[::gpu::platform_other][::gpu::device_other];

      if (::gpu::platform_other == ::gpu::platform &&
          ::gpu::device_other   == ::gpu::device) // only use one device
      {
         vcl::context = cl::Context({vcl::default_device});

         vcl::queue = cl::CommandQueue(vcl::context, vcl::default_device,
                                       CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE);

         vcl::context_other = vcl::context;

         vcl::queue_other = vcl::queue;
      }
      else if (::gpu::platform_other == ::gpu::platform) // use two devices on the same platform
      {
         vcl::context = cl::Context({vcl::default_device, vcl::extra_device});

         vcl::context_other = vcl::context;

         vcl::queue = cl::CommandQueue(vcl::context, vcl::default_device,
                                       CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE);

         vcl::queue_other = cl::CommandQueue(vcl::context, vcl::extra_device,
                                             CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE);
      }
      else // multiple platforms
      {
         vcl::context       = cl::Context({vcl::default_device});
         vcl::context_other = cl::Context({vcl::extra_device});

         vcl::queue = cl::CommandQueue(vcl::context, vcl::default_device,
                                       CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE);

         vcl::queue_other = cl::CommandQueue(vcl::context_other, vcl::extra_device,
                                             CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE);
      }
#endif // ENABLE_MULTIPLE_DEVICES


#ifdef OPENCL_LOG
      vcl::OCLLOG << "Using default platform " << default_platform.getInfo<CL_PLATFORM_NAME>()  << '\n';
      vcl::OCLLOG << "Using default device "   << vcl::default_device.getInfo<CL_DEVICE_NAME>() << '\n';

#ifdef ENABLE_MULTIPLE_DEVICES
      if (::gpu::platform_other != ::gpu::platform ||
          ::gpu::device_other != ::gpu::device)
      {
         vcl::OCLLOG << "Using extra platform " << extra_platform.getInfo<CL_PLATFORM_NAME>()  << '\n';
         vcl::OCLLOG << "Using extra device "   << vcl::extra_device.getInfo<CL_DEVICE_NAME>() << '\n';
      }
#endif // ENABLE_MULTIPLE_DEVICES

      vcl::OCLLOG << '\n';
#endif // OPENCL_LOG

      // build all OpenCL kernels
      vcl::build_kernels();

      success = true;

      success &= vcl::initialize_atoms();
      success &= vcl::initialize_fields();
      success &= vcl::initialize_cells();
      success &= vcl::initialize_materials();
      success &= vcl::initialize_rng();
      success &= vcl::initialize_kernels();

      vcl::queue.finish();

      // initialisation end time
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> diff = end - start;
      std::cout << "OpenCL initialization took " << diff.count() << " seconds." << std::endl;

      // simulation start time
      vcl::time::sim_start = std::chrono::high_resolution_clock::now();
#endif // OPENCL

      return success;

   }

#ifdef OPENCL

   namespace internal
   {
      bool initialize_atoms(void) noexcept
      {
         // Allocate and initialize device memory for atomic spins
         vcl::atoms::spin_array = vcl::Buffer3D(vcl::context, vcl::queue,
                                                CL_MEM_READ_WRITE,
                                                ::atoms::x_spin_array,
                                                ::atoms::y_spin_array,
                                                ::atoms::z_spin_array);

         // Allocate and initialize device memory for atomic coordinates
         vcl::atoms::coord_array = vcl::Buffer3D(vcl::context, vcl::queue,
                                                 CL_MEM_READ_WRITE | CL_MEM_HOST_WRITE_ONLY,
                                                 ::atoms::x_coord_array,
                                                 ::atoms::y_coord_array,
                                                 ::atoms::z_coord_array);

         // Allocate and initialize device memory for atomic information
         vcl::atoms::type_array = vcl::create_device_buffer(::atoms::type_array,
                                                            CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY);

         // Allocate and initialize cell information
         vcl::atoms::cell_array = vcl::create_device_buffer(::atoms::cell_array,
                                                            CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY);

         // Allocate and initialize unrolled spin norm array
         vcl::atoms::spin_norm_array = vcl::create_device_buffer(::atoms::m_spin_array,
                                                                 CL_MEM_READ_WRITE | CL_MEM_HOST_WRITE_ONLY);

         return true;
      }

      bool initialize_fields(void) noexcept
      {
         // Allocate device memory and initialize total spin field arrays
         vcl::total_spin_field_array = vcl::Buffer3D(vcl::context, vcl::queue,
                                                     CL_MEM_READ_WRITE | CL_MEM_HOST_WRITE_ONLY,
                                                     ::atoms::x_total_spin_field_array,
                                                     ::atoms::y_total_spin_field_array,
                                                     ::atoms::z_total_spin_field_array);

         // Allocate device memory and initialize external field arrays
         vcl::total_external_field_array = vcl::Buffer3D(vcl::context, vcl::queue,
                                                         CL_MEM_READ_WRITE | CL_MEM_HOST_WRITE_ONLY,
                                                         ::atoms::x_total_external_field_array,
                                                         ::atoms::y_total_external_field_array,
                                                         ::atoms::z_total_external_field_array);

         // Allocate device memory and initialize for dipolar field
         vcl::dipolar_field_array = vcl::Buffer3D(vcl::context, vcl::queue,
                                                  CL_MEM_READ_WRITE | CL_MEM_HOST_WRITE_ONLY,
                                                  ::atoms::x_dipolar_field_array,
                                                  ::atoms::y_dipolar_field_array,
                                                  ::atoms::z_dipolar_field_array);

         return true;
      }

      bool initialize_cells(void) noexcept
      {
         // Allocate device memory and initialize coordinates
         vcl::cells::coord_array = vcl::Buffer3D(vcl::context, vcl::queue,
                                                 CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY,
                                                 ::cells::x_coord_array,
                                                 ::cells::y_coord_array,
                                                 ::cells::z_coord_array);

         // Allocate device memory and initialize cell magnetization
         vcl::cells::mag_array = vcl::Buffer3D(vcl::context, vcl::queue,
                                               CL_MEM_READ_WRITE | CL_MEM_HOST_WRITE_ONLY,
                                               ::cells::x_mag_array,
                                               ::cells::y_mag_array,
                                               ::cells::z_mag_array);

         // Allocate device memory and initialize cell fields
         vcl::cells::field_array = vcl::Buffer3D(vcl::context, vcl::queue,
                                                 CL_MEM_READ_WRITE | CL_MEM_HOST_WRITE_ONLY,
                                                 ::cells::x_field_array,
                                                 ::cells::y_field_array,
                                                 ::cells::z_field_array);

         // Allocate device memory and initialize voulme array
         vcl::cells::volume_array = vcl::create_device_buffer(::cells::volume_array,
                                                              CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY);

         // Allocate device memory and initialize number of atoms for each cell
         vcl::cells::num_atoms = vcl::create_device_buffer(::cells::num_atoms_in_cell,
                                                           CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY);

         return true;
      }

      bool initialize_materials(void) noexcept
      {
         std::vector<material_parameters_t> h_materials(::mp::num_materials);
         for (unsigned i=0; i<::mp::num_materials; ++i)
         {
            double mu_s_si = ::mp::material[i].mu_s_SI;

            h_materials[i].alpha = ::mp::material[i].alpha;
            h_materials[i].gamma_rel = ::mp::material[i].gamma_rel;
            h_materials[i].mu_s_si = mu_s_si;
            h_materials[i].i_mu_s_si = 1.0 / mu_s_si;
            h_materials[i].k_latt = ::mp::material[i].Klatt_SI / mu_s_si;
            h_materials[i].sh2 = ::mp::material[i].sh2 / mu_s_si;
            h_materials[i].sh4 = ::mp::material[i].sh4 / mu_s_si;
            h_materials[i].sh6 = ::mp::material[i].sh6 / mu_s_si;
            h_materials[i].ku = ::mp::material[i].Ku;
            h_materials[i].anisotropy_unit_x = ::mp::material[i].UniaxialAnisotropyUnitVector[0];
            h_materials[i].anisotropy_unit_y = ::mp::material[i].UniaxialAnisotropyUnitVector[1];
            h_materials[i].anisotropy_unit_z = ::mp::material[i].UniaxialAnisotropyUnitVector[2];
            h_materials[i].applied_field_strength = ::mp::material[i].applied_field_strength;
            h_materials[i].applied_field_unit_x = ::mp::material[i].applied_field_unit_vector[0];
            h_materials[i].applied_field_unit_y = ::mp::material[i].applied_field_unit_vector[1];
            h_materials[i].applied_field_unit_z = ::mp::material[i].applied_field_unit_vector[2];
            h_materials[i].Kc1_SI = ::mp::material[i].Kc1_SI;
            h_materials[i].temperature = ::mp::material[i].temperature;
            h_materials[i].temperature_rescaling_alpha = ::mp::material[i].temperature_rescaling_alpha;
            h_materials[i].temperature_rescaling_Tc = ::mp::material[i].temperature_rescaling_Tc;
            h_materials[i].H_th_sigma = ::mp::material[i].H_th_sigma;
         }

         // Allocate device memory and initialize materials array
         vcl::mp::materials = vcl::create_device_buffer(h_materials,
                                                        CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY,
                                                        CL_TRUE);

         return true;
      }

      bool initialize_rng(void) noexcept
      {
         // each atom needs three random numbers per Heun step
         // Box Muller transform uses two uniform random numbers to
         // generate two gaussian random numbers so an even number of
         // random numbers need to be stored
         if (::atoms::num_atoms % 2 == 0)
         {
            vcl::rng::n_rands = ::atoms::num_atoms*3;
         }
         else
         {
            vcl::rng::n_rands = ::atoms::num_atoms*3 + 1;
         }

         std::vector<cl_ulong> rs(vcl::rng::n_rands);

         const size_t g_buffer_size = vcl::rng::n_rands * sizeof (vcl::real_t);

         std::random_device rd;
         const std::mt19937_64 gen(rd());

         // must not seed xorshift with 0
         const std::uniform_int_distribution<cl_ulong> dist(1);

         auto rand64 = std::bind(dist, gen);

         for (unsigned i=0; i<vcl::rng::n_rands; ++i)
         {
            rs[i] = rand64();
         }

#ifdef ENABLE_MULTIPLE_DEVICES
         vcl::rng::state = vcl::create_device_buffer(rs,
                                                     CL_MEM_READ_WRITE | CL_MEM_HOST_WRITE_ONLY,
                                                     CL_TRUE,
                                                     vcl::queue_other,
                                                     vcl::context_other);

         if (::gpu::platform_other != ::gpu::platform)
         {
            // need to read grands on the host
            vcl::rng::grands = cl::Buffer(vcl::context_other,
                                          CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY,
                                          g_buffer_size);

            // copy of grands on main platform
            vcl::rng::grands_copy = cl::Buffer(vcl::context,
                                               CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY,
                                               g_buffer_size);
         }
         else
         {
            vcl::rng::grands = cl::Buffer(vcl::context_other,
                                          CL_MEM_READ_WRITE,
                                          g_buffer_size);
         }
                                          
#else
         vcl::rng::grands = cl::Buffer(vcl::context,
                                       CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS,
                                       g_buffer_size);

         vcl::rng::state = vcl::create_device_buffer(rs,
                                                     CL_MEM_READ_WRITE | CL_MEM_HOST_WRITE_ONLY,
                                                     CL_TRUE,
                                                     vcl::queue,
                                                     vcl::context);
#endif // ENABLE_MULTIPLE_DEVICES

         return true;
      }
   }

#endif // OPENCL

} // end of vopencl namespace

