//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>

#ifdef OPENCL_DEBUG
#include <chrono>
#endif // OPENCL_DEBUG

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

#ifdef ENABLE_OPENCL_TESTS
#include "tests/tests.hpp"
#endif

#ifdef OPENCL
namespace vcl = ::vopencl::internal;
#endif

namespace vopencl
{
   namespace internal
   {
      cl::NDRange global;
   }

   //----------------------------------------------------------------------------
   // Function to initialize vopencl module
   //----------------------------------------------------------------------------
   bool initialize(bool cpu_stats)
   {
      bool success = false;

#ifdef OPENCL

#ifdef OPENCL_DEBUG
      auto start = std::chrono::high_resolution_clock::now();
#endif //OPENCL_DEBUG

      std::string message("OpenCL has been enabled in ");
#ifdef OPENCL_DP
      message.append("double precision mode.");
#else
      message.append("single precision mode.");
#endif // OPENCL_DP

#ifdef OPENCL_USE_NATIVE_FUNCTIONS
      message.append(" Native functions will be used.");
#endif

      std::cout << message << std::endl;
      zlog << zTs() << message << std::endl;

      //vcl::stats::use_cpu = cpu_stats;
      vcl::stats::use_cpu = true;

      // find OpenCL platforms and devices
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

      std::vector<std::vector<cl::Device>> devices(nplatforms);
      unsigned ndevices = 0;
      for (unsigned i=0; i<nplatforms; ++i)
      {
         std::vector<cl::Device> tmp_devices;
         platforms[i].getDevices(CL_DEVICE_TYPE_ALL, &tmp_devices);
         devices[i] = tmp_devices;
         ndevices += tmp_devices.size();

#ifdef OPENCL_DEBUG
         vcl::OCLLOG << "Found platform " << platforms[i].getInfo<CL_PLATFORM_NAME>() << std::endl;
         for (unsigned j=0; j<tmp_devices.size(); ++j)
         {
            vcl::OCLLOG << "Found device " << tmp_devices[j].getInfo<CL_DEVICE_NAME>() << std::endl;
            vcl::OCLLOG << "with version " << tmp_devices[j].getInfo<CL_DEVICE_VERSION>() << std::endl;
         }
#endif // OPENCL_DEBUG
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

      cl::Platform default_platform = platforms[0];
      vcl::default_device = devices[0][0];

#ifdef OPENCL_DEBUG
      vcl::OCLLOG << "Using default platform " << default_platform.getInfo<CL_PLATFORM_NAME>() << std::endl;
      vcl::OCLLOG << "Using default device " << vcl::default_device.getInfo<CL_DEVICE_NAME>() << std::endl;
#endif // OPENCL_DEBUG

      vcl::context = cl::Context({vcl::default_device});

      vcl::queue = cl::CommandQueue(vcl::context, vcl::default_device, CL_QUEUE_PROFILING_ENABLE);

      if (::gpu::num_threads > 0)
      {
         vcl::global = cl::NDRange(::gpu::num_threads);
      }
      else
      {
         int default_nthreads = 64;
         terminaltextcolor(YELLOW);
         std::cerr << "Warning: defaulting to " << default_nthreads << " threads." << std::endl;
         std::cerr << "Use gpu:num-threads=n in the input file to specify number." << std::endl;
         terminaltextcolor(WHITE);
         vcl::global = cl::NDRange(default_nthreads);
      }

      // build all OpenCL kernels
      vcl::build_kernels();

#ifdef ENABLE_OPENCL_TESTS
      if(!vcl::test::all())
      {
         terminaltextcolor(RED);
         std::cerr << "At least one OpenCL test has failed. Aborting." << std::endl;
         terminaltextcolor(WHITE);
         ::err::vexit();
      }
#endif

      success = true;

      success &= vcl::initialize_atoms();
      success &= vcl::initialize_fields();
      success &= vcl::initialize_cells();
      success &= vcl::initialize_materials();
      success &= vcl::initialize_topology();
      //success &= vcl::initialize_stats();
      success &= vcl::initialize_rng();
      success &= vcl::initialize_kernels();

      vcl::queue.finish();

#ifdef OPENCL_DEBUG
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> diff = end - start;
      std::cout << "OpenCL initialization took " << diff.count() << " seconds." << std::endl;
#endif // OPENCL_DEBUG
#endif // OPENCL

      return success;

   }

#ifdef OPENCL

   namespace internal
   {
      bool initialize_atoms(void) noexcept
      {
         // Allocate and initialize device memory for atomic spins
         vcl::atoms::spin_array = vcl::Buffer3D<vcl::real_t>(vcl::context, vcl::queue,
                                                             CL_MEM_READ_WRITE,
                                                             ::atoms::x_spin_array,
                                                             ::atoms::y_spin_array,
                                                             ::atoms::z_spin_array);

         // Allocate and initialize device memory for atomic coordinates
         vcl::atoms::coord_array = vcl::Buffer3D<vcl::real_t>(vcl::context, vcl::queue,
                                                              CL_MEM_READ_WRITE,
                                                              ::atoms::x_coord_array,
                                                              ::atoms::y_coord_array,
                                                              ::atoms::z_coord_array);

         // Allocate and initialize device memory for atomic information
         size_t buff_size = ::atoms::type_array.size() * sizeof ::atoms::type_array[0];
         vcl::atoms::type_array = cl::Buffer(vcl::context, CL_MEM_READ_ONLY, buff_size);
         vcl::queue.enqueueWriteBuffer(vcl::atoms::type_array,
                                       CL_FALSE,
                                       0,
                                       buff_size,
                                       ::atoms::type_array.data());

         // Allocate and initialize cell information
         buff_size = ::atoms::cell_array.size() * sizeof ::atoms::cell_array[0];
         vcl::atoms::cell_array = cl::Buffer(vcl::context, CL_MEM_READ_ONLY, buff_size);
         vcl::queue.enqueueWriteBuffer(vcl::atoms::cell_array,
                                       CL_FALSE,
                                       0,
                                       buff_size,
                                       ::atoms::cell_array.data());

         // Allocate and initialize unrolled spin norm array
         buff_size = ::atoms::m_spin_array.size() * sizeof ::atoms::m_spin_array[0];
         vcl::atoms::spin_norm_array = cl::Buffer(vcl::context, CL_MEM_READ_WRITE, buff_size);
         vcl::queue.enqueueWriteBuffer(vcl::atoms::spin_norm_array,
                                       CL_FALSE,
                                       0,
                                       buff_size,
                                       ::atoms::m_spin_array.data());

         return true;
      }

      bool initialize_fields(void) noexcept
      {
         // Allocate device memory and initialize total spin field arrays
         vcl::total_spin_field_array = vcl::Buffer3D<vcl::real_t>(vcl::context, vcl::queue, CL_MEM_READ_WRITE,
                                                                  ::atoms::x_total_spin_field_array,
                                                                  ::atoms::y_total_spin_field_array,
                                                                  ::atoms::z_total_spin_field_array);

         // Allocate device memory and initialize external field arrays
         vcl::total_external_field_array = vcl::Buffer3D<vcl::real_t>(vcl::context, vcl::queue, CL_MEM_READ_WRITE,
                                                                      ::atoms::x_total_external_field_array,
                                                                      ::atoms::y_total_external_field_array,
                                                                      ::atoms::z_total_external_field_array);

         // Allocate device memory and initialize for dipolar field
         vcl::dipolar_field_array = vcl::Buffer3D<vcl::real_t>(vcl::context, vcl::queue, CL_MEM_READ_WRITE,
                                                               ::atoms::x_dipolar_field_array,
                                                               ::atoms::y_dipolar_field_array,
                                                               ::atoms::z_dipolar_field_array);

         return true;
      }

      bool initialize_cells(void) noexcept
      {
         // Allocate device memory and initialize coordinates
         vcl::cells::coord_array = vcl::Buffer3D<vcl::real_t>(vcl::context, vcl::queue, CL_MEM_READ_ONLY,
                                                              ::cells::x_coord_array,
                                                              ::cells::y_coord_array,
                                                              ::cells::z_coord_array);

         // Allocate device memory and initialize cell magnetization
         vcl::cells::mag_array = vcl::Buffer3D<vcl::real_t>(vcl::context, vcl::queue, CL_MEM_READ_WRITE,
                                                            ::cells::x_mag_array,
                                                            ::cells::y_mag_array,
                                                            ::cells::z_mag_array);

         // Allocate device memory and initialize cell fields
         vcl::cells::field_array = vcl::Buffer3D<vcl::real_t>(vcl::context, vcl::queue, CL_MEM_READ_WRITE,
                                                              ::cells::x_field_array,
                                                              ::cells::y_field_array,
                                                              ::cells::z_field_array);

         // Allocate device memory and initialize voulme array
         size_t buff_size = ::cells::volume_array.size() * sizeof ::cells::volume_array[0];
         vcl::cells::volume_array = cl::Buffer(vcl::context, CL_MEM_READ_ONLY, buff_size);
         vcl::queue.enqueueWriteBuffer(vcl::cells::volume_array,
                                       CL_FALSE,
                                       0,
                                       buff_size,
                                       ::cells::volume_array.data());

         // Allocate device memory and initialize number of atoms for each cell
         buff_size = ::cells::num_atoms_in_cell.size() * sizeof ::cells::num_atoms_in_cell[0];
         vcl::cells::num_atoms = cl::Buffer(vcl::context, CL_MEM_READ_ONLY, buff_size);
         vcl::queue.enqueueWriteBuffer(vcl::cells::num_atoms,
                                       CL_FALSE,
                                       0,
                                       buff_size,
                                       ::cells::num_atoms_in_cell.data());

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
         const size_t buff_size = ::mp::num_materials * sizeof h_materials[0];
         vcl::mp::materials = cl::Buffer(vcl::context, CL_MEM_READ_ONLY, buff_size);
         vcl::queue.enqueueWriteBuffer(vcl::mp::materials,
                                       CL_FALSE,
                                       0,
                                       buff_size,
                                       h_materials.data());

         return true;
      }

      bool initialize_topology(void) noexcept
      {
         std::vector<cl_uint> limits_h(::atoms::num_atoms+1);
         limits_h[0] = 0;
         for (int atom=0; atom<::atoms::num_atoms; ++atom)
         {
            limits_h[atom+1] = ::atoms::neighbour_list_end_index[atom]+1;
         }

         // Allocate device memory and initialize limits array
         size_t buff_size = limits_h.size() * sizeof limits_h[0];
         vcl::atoms::limits = cl::Buffer(vcl::context, CL_MEM_READ_ONLY, buff_size);
         vcl::queue.enqueueWriteBuffer(vcl::atoms::limits,
                                       CL_FALSE,
                                       0,
                                       buff_size,
                                       limits_h.data());

         buff_size = ::atoms::neighbour_list_array.size() * sizeof ::atoms::neighbour_list_array[0];
         vcl::atoms::neighbours = cl::Buffer(vcl::context, CL_MEM_READ_ONLY, buff_size);
         vcl::queue.enqueueWriteBuffer(vcl::atoms::neighbours,
                                       CL_FALSE,
                                       0,
                                       buff_size,
                                       ::atoms::neighbour_list_array.data());

         return true;
      }


      static cl_ulong rand64(void) noexcept
      {
         cl_ulong r = std::rand();
         return (r << 32) | std::rand();
      }

      bool initialize_rng(void) noexcept
      {
         // each atom needs three random numbers per Heun step
         // Box Muller transform uses two uniform random numbers to
         // generate two gaussian random numbers so an even number of
         // random numbers need to be stored
         size_t n_rands;
         if (::atoms::num_atoms % 2 == 0)
            n_rands = ::atoms::num_atoms*3;
         else
            n_rands = ::atoms::num_atoms*3 + 1;

         std::vector<cl_ulong> rs(n_rands);

         const size_t u_buffer_size = rs.size() * sizeof(cl_ulong);
         const size_t g_buffer_size = rs.size() * sizeof(vcl::real_t);

         vcl::rng::state  = cl::Buffer(vcl::context, CL_MEM_READ_WRITE, u_buffer_size);
         vcl::rng::grands = cl::Buffer(vcl::context, CL_MEM_READ_WRITE, g_buffer_size);

         std::srand(std::time(NULL));
         for (unsigned i=0; i<rs.size(); ++i)
         {
            // must not seed xorshift with 0
            cl_ulong r;
            do { r = rand64(); } while (r == 0);
            assert(r != cl_ulong(0));
            rs[i] = r;
         }

         vcl::queue.enqueueWriteBuffer(vcl::rng::state, CL_FALSE, 0, u_buffer_size, rs.data());

         return true;
      }
   }

#endif // OPENCL

} // end of vopencl namespace

