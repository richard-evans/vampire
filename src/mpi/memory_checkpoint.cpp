//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2026. All rights reserved.
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <fstream>
#include <string>

// System headers for the current resident set size
#if defined(__APPLE__)
   #include <mach/mach.h>
#elif defined(__linux__)
   #include <unistd.h>
#endif

// Vampire headers
#include "vio.hpp"
#include "vmpi.hpp"
#include "vutil.hpp"

namespace vmpi{

namespace{

//------------------------------------------------------------------------------
// Current (not peak) resident set size of this process in kB: /proc/self/statm on Linux, the task's resident size
// on macOS, -1 elsewhere
//------------------------------------------------------------------------------
double current_memory_usage_kb(){

   #if defined(__APPLE__)
      mach_task_basic_info_data_t info;
      mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
      if(task_info(mach_task_self(), MACH_TASK_BASIC_INFO, reinterpret_cast<task_info_t>(&info), &count) == KERN_SUCCESS){
         return double(info.resident_size) / 1024.0;
      }
      return -1.0;
   #elif defined(__linux__)
      std::ifstream statm("/proc/self/statm");
      long pages_total = 0, pages_resident = 0;
      if(statm >> pages_total >> pages_resident) return double(pages_resident) * double(sysconf(_SC_PAGESIZE)) / 1024.0;
      return -1.0;
   #else
      return -1.0;
   #endif

}

}

//------------------------------------------------------------------------------
// Collective: writes one MEMORY_CHECKPOINT line to the log with the peak and current resident set size (and the
// proportional set size on Linux, which counts MPI shared memory once per region) of the worst rank and the mean
// over ranks, in MB. Peak RSS only ever rises, so the checkpoint at which the worst-rank peak jumps locates the
// set-up phase responsible, and a rise in peak without a matching rise in current RSS shows a transient.
//------------------------------------------------------------------------------
void memory_checkpoint(const std::string& label){

   vmpi::barrier();

   const double peak_mb    = vutil::peak_memory_usage_kb() / 1024.0;
   const double current_mb = current_memory_usage_kb() / 1024.0;
   const double pss_kb     = vutil::proportional_memory_usage_kb();
   const double pss_mb     = pss_kb < 0.0 ? -1.0 : pss_kb / 1024.0; // -1 where unavailable (not Linux)

   const double peak_max    = vmpi::reduce_max(peak_mb);
   const double current_max = vmpi::reduce_max(current_mb);
   const double pss_max     = vmpi::reduce_max(pss_mb);
   const double nprocs      = double(vmpi::num_processors > 0 ? vmpi::num_processors : 1);
   const double peak_mean    = vmpi::reduce_sum(peak_mb) / nprocs;
   const double current_mean = vmpi::reduce_sum(current_mb) / nprocs;
   const double pss_mean     = vmpi::reduce_sum(pss_mb) / nprocs;

   if(vmpi::master){
      zlog << zTs() << "MEMORY_CHECKPOINT label=" << label
           << " peak_max_mb=" << peak_max << " peak_mean_mb=" << peak_mean
           << " current_max_mb=" << current_max << " current_mean_mb=" << current_mean
           << " pss_max_mb=" << pss_max << " pss_mean_mb=" << pss_mean << std::endl;
   }

}

}
