#!/bin/bash
#------------------------------------------------------------------------------
#
#   This file is part of the VAMPIRE open source package under the
#   Free BSD licence (see licence file for details).
#
#------------------------------------------------------------------------------
# Builds the CUDA version of vampire in double and single precision, runs the
# validation suite with each build and compares against the CPU references.
# Each build is also run with statistics calculated on the CPU
# (gpu:calculate-statistics-on-cpu=true, two seeds per statistical case): the
# dynamics are the same, so every run must agree with the same run of the
# default build to rounding of the statistics sums. Also checks that the CPU
# code path of the CUDA binary (gpu:acceleration=off, one seed per case)
# reproduces the serial references byte for byte, and times the hamr-perf
# benchmark. Finally runs the HAMR chain of hamr_chain.sh: the CPU path relaxes
# the medium by Monte Carlo and writes it as a reference, each build writes it
# from the same relaxed state (compared with the reference within the
# deterministic tolerance), and the CPU path loads the checkpoint of each GPU
# write, whose spins must be bit for bit those of that write.
#
# Usage (from anywhere):  developers/test/validation/run_gpu_suite.sh [label-prefix]
#
# Environment variables:
#   GPU_JOBS   validation runs sharing the GPU at the same time (default: 4);
#              results do not depend on it, benchmarks always run alone
#
# Results are written to developers/test/validation/results/<prefix>-* and packed, with
# build logs, comparison reports and system information, into
# developers/test/validation/results/gpu-suite-<host>-<date>.tar.gz
#
# Note: the CUDA and serial builds share object files, so this script runs
# 'make clean' before each build and again at the end.
#------------------------------------------------------------------------------
set -u
cd "$(dirname "$0")/../../.."

prefix=${1:-gpu}
host=$(hostname -s)
name=gpu-suite-$host-$(date +%Y%m%d-%H%M)
results=developers/test/validation/results
out=$results/$name
vt="python3 developers/test/validation/vtest.py"
source developers/test/validation/hamr_chain.sh
gpu_jobs=${GPU_JOBS:-4}
# CPU code path: every case except benchmarks and the granular cases, whose (long) CPU runs
# test no code that differs between the builds
cpu_cases=$(cd developers/test/validation/cases && for c in */; do c=${c%/}; grep -q '"benchmark"' $c/case.json || echo $c; done | grep -v '^hamr-granular' | paste -sd,)
mkdir -p $out

# record system information
{ nvidia-smi; nvcc --version; g++ --version; echo "git: $(git rev-parse HEAD)"; git status --short; } > $out/system.txt 2>&1

labels=""
for pr in dp sp; do

   echo "=== Building CUDA version ($pr) ==="
   make clean
   if ! make -j$(nproc) cuda CUDA_PR=-DCUDA_$(echo $pr | tr a-z A-Z) > $out/build-$pr.log 2>&1; then
      echo "Build failed, see $out/build-$pr.log"
      continue
   fi
   cp vampire-cuda vampire-cuda-$pr

   echo "=== Running validation suite ($pr) ==="
   $vt run --exe ./vampire-cuda-$pr --label $prefix-$pr --precision $pr --jobs $gpu_jobs
   $vt compare --test $prefix-$pr --ref cpu --verbose > $out/compare-$pr.txt

   echo "=== Running validation suite ($pr, statistics on CPU) ==="
   $vt run --exe ./vampire-cuda-$pr --label $prefix-$pr-cpustats --precision $pr --seeds 2 --jobs $gpu_jobs --extra "gpu:calculate-statistics-on-cpu=true"
   # same dynamics as the default run, only the statistics summation order differs
   $vt compare --test $prefix-$pr-cpustats --ref $prefix-$pr --paired 1e-10 --verbose > $out/compare-$pr-cpustats.txt

   echo "=== Running hamr-perf benchmark ($pr) ==="
   $vt run --exe ./vampire-cuda-$pr --label $prefix-$pr-bench --precision $pr --cases hamr-perf
   $vt compare --test $prefix-$pr-bench --ref cpu --cases hamr-perf > $out/benchmark-$pr.txt
   $vt run --exe ./vampire-cuda-$pr --label $prefix-$pr-cpustats-bench --precision $pr --cases hamr-perf --extra "gpu:calculate-statistics-on-cpu=true"
   $vt compare --test $prefix-$pr-cpustats-bench --ref cpu --cases hamr-perf > $out/benchmark-$pr-cpustats.txt

   labels="$labels $prefix-$pr $prefix-$pr-bench $prefix-$pr-cpustats $prefix-$pr-cpustats-bench"
done

# the CPU path of the CUDA binary should be bit-identical to vampire-serial
if [ -x vampire-cuda-dp ]; then
   echo "=== Running validation suite (CUDA binary, gpu:acceleration=off) ==="
   $vt run --exe ./vampire-cuda-dp --label $prefix-cpupath --cases $cpu_cases --seeds 1 --extra "gpu:acceleration=off" --jobs $(nproc)
   $vt identical --test $prefix-cpupath --ref cpu > $out/compare-cpupath.txt
   labels="$labels $prefix-cpupath"

   echo "=== HAMR chain: CPU relaxation and reference write, GPU writes, CPU loads of their checkpoints ==="
   cpu="gpu:acceleration=off"
   chain_run relax - $prefix-chain-relax dp "$cpu" $PWD/vampire-cuda-dp
   chain_run write $prefix-chain-relax $prefix-chain-cpu dp "$cpu" $PWD/vampire-cuda-dp
   labels="$labels $prefix-chain-relax $prefix-chain-cpu"
   for pr in dp sp; do
      [ -x vampire-cuda-$pr ] || continue
      chain_run write $prefix-chain-relax $prefix-$pr-chain $pr "" $PWD/vampire-cuda-$pr
      chain_run load $prefix-$pr-chain $prefix-$pr-chainload $pr "$cpu" $PWD/vampire-cuda-dp
      {
         $vt compare --test $prefix-$pr-chain --ref $prefix-chain-cpu --cases hamr-granular-t0 --verbose
         $vt configs --test $prefix-$pr-chainload --ref $prefix-$pr-chain --cases hamr-granular-t0
      } > $out/chain-$pr.txt
      labels="$labels $prefix-$pr-chain $prefix-$pr-chainload"
   done

   cp vampire-cuda-dp vampire-cuda   # leave the default (double precision) build in place
fi
make clean

# configuration files must never contain NaN: on the GPU the spins on the host are
# overwritten with NaN whenever they are not valid, so NaN here means invalid spins were used
grep -rli nan --include='*.data' --include='*.cfg' $(for l in $labels; do echo $results/$l; done) > $out/config-nan-check.txt 2>/dev/null
echo "configuration files containing NaN: $(wc -l < $out/config-nan-check.txt)" >> $out/config-nan-check.txt

# pack results (atomic configuration files are excluded to keep the archive small)
tar czf $results/$name.tar.gz --exclude='*.data' --exclude='*.meta' --exclude='*.chk' -C $results $name $labels

echo
echo "=== Summary ==="
for f in $out/compare-*.txt $out/benchmark-*.txt $out/config-nan-check.txt; do
   [ -f $f ] && echo "$(basename $f): $(tail -n 1 $f)"
done
for f in $out/chain-*.txt; do
   echo "$(basename $f) (GPU write vs CPU write; CPU load vs GPU write): $(grep OVERALL $f | paste -sd';')"
done
echo "Archive: $results/$name.tar.gz"
