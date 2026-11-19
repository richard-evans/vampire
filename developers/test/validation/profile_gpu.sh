#!/bin/bash
#------------------------------------------------------------------------------
#
#   This file is part of the VAMPIRE open source package under the
#   Free BSD licence (see licence file for details).
#
#------------------------------------------------------------------------------
# Profiles the CUDA version of vampire with Nsight Systems on one validation
# case (default: hamr-perf).
#
#   1. builds vampire-cuda (with -lineinfo) in the requested precision
#   2. runs the case once without the profiler, logging wall time and GPU
#      utilisation every second (nvidia-smi dmon)
#   3. runs the case under nsys, capturing a time window after initialisation
#      (CUDA API/kernels/memory copies, OS runtime calls and CPU sampling);
#      the run is stopped when the window ends
#   4. exports summary tables (csv) and an sqlite database of the window
#   5. optionally (NCU=1) profiles a few launches of the main kernels with
#      Nsight Compute, which shows whether each kernel is limited by memory
#      bandwidth, memory latency or arithmetic
#
# Usage (from anywhere):  developers/test/validation/profile_gpu.sh [case] [dp|sp]
#
# Environment variables:
#   NSYS_DELAY     seconds to wait before capturing (default 30, skips set-up)
#   NSYS_DURATION  length of the captured window in seconds (default 20)
#   JOBS           parallel build jobs (default: number of cores)
#   EXTRA          line appended to the input file, for example
#                  EXTRA="gpu:calculate-statistics-on-cpu=true"
#   NCU            set to 1 to add the Nsight Compute step
#   NCU_KERNELS    kernel name regular expression (default: the field, Heun,
#                  thermal and statistics kernels)
#   NCU_SKIP       matching launches to skip first (default 15000, past the
#                  hamr-perf equilibration)
#   NCU_COUNT      matching launches to profile (default 12)
#
# Nsight Compute needs access to the GPU performance counters; if it reports
# ERR_NVGPUCTRPERM, run the script with sudo or enable counter access for
# all users in the NVIDIA driver.
#
# Results are packed into developers/test/validation/results/profile-<case>-<pr>-<host>-<date>.tar.gz
#
# Note: the CUDA and serial builds share object files, so this script runs
# 'make clean' before and after the build.
#------------------------------------------------------------------------------
set -u
cd "$(dirname "$0")/../../.."

case=${1:-hamr-perf}
pr=${2:-dp}
delay=${NSYS_DELAY:-30}
duration=${NSYS_DURATION:-20}
jobs=${JOBS:-$(nproc)}

name=profile-$case-$pr-$(hostname -s)-$(date +%Y%m%d-%H%M)
results=developers/test/validation/results
out=$results/$name
abs_out=$PWD/$out
mkdir -p $out
exe=$PWD/vampire-cuda-$pr-profile

command -v nsys > /dev/null || { echo "nsys not found in PATH"; exit 1; }
[ -d developers/test/validation/cases/$case ] || { echo "unknown case: $case"; exit 1; }

# record system information
{ echo "extra input: ${EXTRA:-none}"; nvidia-smi; nvcc --version; nsys --version; g++ --version; echo "git: $(git rev-parse HEAD)"; git status --short; } > $out/system.txt 2>&1

echo "=== Building CUDA version ($pr, $jobs jobs) ==="
make clean > /dev/null
# -lineinfo maps kernel time to source lines; it does not change the generated code
if ! make -j$jobs cuda CUDA_PR="-DCUDA_$(echo $pr | tr a-z A-Z) -lineinfo" > $out/build.log 2>&1; then
   echo "Build failed, see $out/build.log"
   exit 1
fi
mv vampire-cuda $exe
make clean > /dev/null

# each run gets a fresh copy of the case directory
run_dir(){
   rm -rf $out/$1; cp -r developers/test/validation/cases/$case $out/$1; rm -f $out/$1/case.json
   if [ -n "${EXTRA:-}" ]; then echo "$EXTRA" >> $out/$1/input; fi
}

echo "=== Unprofiled run (wall time and GPU utilisation) ==="
run_dir plain
nvidia-smi dmon -s um -d 1 -o T > $out/dmon.txt 2>&1 &
dmon=$!
start=$(date +%s.%N)
( cd $out/plain && $exe > log.txt 2>&1 )
end=$(date +%s.%N)
kill $dmon 2> /dev/null
awk -v s=$start -v e=$end 'BEGIN { printf "wall time: %.1f s\n", e - s }' | tee $out/walltime.txt

echo "=== Profiled run (window: ${delay}s + ${duration}s) ==="
run_dir nsys
( cd $out/nsys && nsys profile --trace=cuda,osrt --sample=process-tree \
     --delay=$delay --duration=$duration \
     --output=$abs_out/report --force-overwrite=true $exe > log.txt 2>&1 )

echo "=== Exporting reports ==="
nsys stats --report cuda_gpu_kern_sum,cuda_gpu_mem_time_sum,cuda_gpu_mem_size_sum,cuda_api_sum,osrt_sum,cuda_gpu_trace \
   --format csv --output $out/stats $out/report.nsys-rep > $out/stats.log 2>&1
nsys export --type sqlite --output $out/report.sqlite --force-overwrite=true $out/report.nsys-rep >> $out/stats.log 2>&1

if [ -n "${NCU:-}" ]; then
   echo "=== Nsight Compute (kernels: ${NCU_KERNELS:-default}) ==="
   command -v ncu > /dev/null || { echo "ncu not found in PATH"; exit 1; }
   run_dir ncu
   ( cd $out/ncu && ncu --kernel-name "regex:${NCU_KERNELS:-update_isotropic|apply_local_temperature|apply_global_temperature|llg_heun|sum_magnetization}" \
        --launch-skip ${NCU_SKIP:-15000} --launch-count ${NCU_COUNT:-12} --set detailed --kill yes \
        -o $abs_out/kernels -f $exe > log.txt 2>&1 )
   ncu --import $(ls $out/kernels.ncu-rep*) --page details --csv > $out/ncu_details.csv 2> $out/ncu_export.log
   tail -n 3 $out/ncu/log.txt
   # ncu reports failures (e.g. no counter permission) in the run log, not via its exit code
   grep -h '==ERROR==' $out/ncu/log.txt $out/ncu_export.log && echo "WARNING: Nsight Compute step failed (see above)"
fi

# pack results (atomic configuration files are excluded to keep the archive small)
tar czf $results/$name.tar.gz --exclude='*.data' --exclude='*.meta' -C $results $name
rm -f $exe

echo
echo "=== Summary ==="
cat $out/walltime.txt
awk '!/^#/ && $3 ~ /^[0-9]+$/ { s += $3; n++ } END { if (n) printf "mean SM utilisation (unprofiled run): %.0f %%\n", s/n }' $out/dmon.txt
echo "Archive: $results/$name.tar.gz ($(du -h $results/$name.tar.gz | cut -f1))"
