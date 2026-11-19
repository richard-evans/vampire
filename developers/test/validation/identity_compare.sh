#!/bin/bash
#------------------------------------------------------------------------------
#
#   This file is part of the VAMPIRE open source package under the
#   Free BSD licence (see licence file for details).
#
#------------------------------------------------------------------------------
# Checks that the CUDA version built from the working tree gives results
# identical, byte for byte, to the CUDA version built from a reference commit,
# for changes that are not meant to change any result.
#
#   1. builds both versions in double and single precision
#   2. runs every validation case with both versions, once as committed and
#      once with extra input lines that exercise the statistics output:
#      output every 150 steps, time averaged, susceptibility, standard
#      deviation, height and grain statistics, configuration files and
#      checkpoints during the simulation
#   3. compares every file written by each run (output, grain, configuration
#      and checkpoint files); comment lines starting with '#' (date, process
#      id, path, git hash) and the stdout and log files (timings) are ignored
#   4. runs the cases with non-magnetic (keep) atoms with and without
#      gpu:skip-non-magnetic-atoms using the working tree version: output and
#      grain files must be identical, and so must the final spins of every
#      magnetic atom
#   5. times the hamr-perf benchmark and a small (5 nm) thermal-perf system
#
# With MGPU_RANKS=N the multi-GPU version (make cuda-multigpu) is compared
# instead, every run on N MPI ranks (which may share one GPU), and the timed
# benchmark is the granular HAMR medium (32 x 32 x 20 nm, periodic in x and y,
# 2000 time steps at 300 K) on one rank.
#
# Usage (from anywhere):  developers/test/validation/identity_compare.sh [ref-commit]
#                         (default: develop)
#
# Environment variables:
#   SEEDS     seeds per statistical case (default: 1; a byte for byte comparison needs
#             only one, statistical comparisons are for run_gpu_suite.sh and ab_compare.sh)
#   JOBS      parallel build jobs (default: number of cores)
#   GPU_JOBS  validation runs sharing the GPU at the same time (default: 4); results do
#             not depend on it, the timed benchmarks always run alone
#   MGPU_RANKS compare the multi-GPU version on this number of ranks (see above)
#   MPIRUN    MPI launcher for MGPU_RANKS (default: mpirun --oversubscribe --bind-to none)
#   NO_CACHE  set to 1 to rebuild and rerun the reference version even if it is cached
#   TIMINGS   set to 0 to skip the timing of the reference version: the benchmark cases
#             still run with the working tree version and their outputs are still
#             compared, against the cached reference outputs (for changes not meant to
#             affect speed; the default 1 times both versions under the same conditions)
#
# Results are packed into developers/test/validation/results/identity-<ref>-<host>-<date>.tar.gz
#
# The reference binaries and the outputs of their validation runs are kept in
# results/identity-cache and reused by later runs with the same reference commit,
# computer, GPU and driver, CUDA compiler, build settings, number of ranks, seeds,
# validation cases, vtest.py and extra input lines (the two most recent sets are
# kept). The working tree version is always built and run in full; the timed
# benchmarks run both versions again (unless TIMINGS=0).
#
# Note: the reference version is built from an exported copy of the commit; the
# working tree version is built incrementally in results/identity-build (see
# incremental_build.sh); the working tree itself is not cleaned.
#------------------------------------------------------------------------------
set -u
cd "$(dirname "$0")/../../.."

ref=${1:-develop}
seeds=${SEEDS:-1}
jobs=${JOBS:-$(nproc)}
gpu_jobs=${GPU_JOBS:-4}
sha=$(git rev-parse --short "$ref") || exit 1
name=identity-$sha-$(hostname -s)-$(date +%Y%m%d-%H%M)
results=developers/test/validation/results
out=$results/$name
mkdir -p $out
vt="python3 developers/test/validation/vtest.py"

# version to compare: the CUDA version, or the multi-GPU version on MGPU_RANKS ranks
target=cuda; binary=vampire-cuda; ranks=()
if [ -n "${MGPU_RANKS:-}" ]; then
   target=cuda-multigpu; binary=vampire-cuda-multigpu; ranks=(--ranks $MGPU_RANKS)
   export MPIRUN=${MPIRUN:-mpirun --oversubscribe --bind-to none}
fi

# extra input lines for the second run of each case
stress=(
   "output:output-rate=150"
   "output:mean-magnetisation"
   "output:mean-magnetisation-length"
   "output:material-mean-magnetisation"
   "output:mean-height-magnetisation"
   "output:mean-susceptibility"
   "output:material-mean-susceptibility"
   "output:material-standard-deviation"
   "grain:output-rate=150"
   "grain:mean-magnetisation-length"
   "grain:mean-susceptibility"
   "config:atoms"
   "config:atoms-output-rate=211"   # every 211 calls of the output routine (every step in HAMR)
   "sim:save-checkpoint=continuous"
   "sim:save-checkpoint-rate=1000"
)
extra=(); for line in "${stress[@]}"; do extra+=(--extra "$line"); done

# cases with atoms of non-magnetic=keep materials, run with and without gpu:skip-non-magnetic-atoms
keep_cases=hamr-granular,hamr-granular-t0
configs=(--extra "config:atoms=end" --extra "config:output-format=binary")

# timed benchmarks: hamr-perf and a 5 nm thermal-perf system (host overhead bound)
small=(--extra "dimensions:system-size-x = 5.0 !nm" --extra "dimensions:system-size-y = 5.0 !nm")
benches=("bench hamr-perf" "bench-small thermal-perf")
if [ -n "${MGPU_RANKS:-}" ]; then
   small=(--ranks 1 --seeds 1 --extra "dimensions:system-size-x = 32 !nm" --extra "dimensions:system-size-y = 32 !nm"
          --extra "sim:program = time-series" --extra "sim:temperature = 300.0" --extra "sim:total-time-steps = 2000")
   benches=("bench-small hamr-granular")
fi

#------------------------------------------------------------------------------
# cache of the reference version: everything its outputs depend on is in the key
#------------------------------------------------------------------------------
key_text=$(
   echo "ref $(git rev-parse "$ref") host $(hostname -s) target $target ranks ${MGPU_RANKS:-} seeds $seeds"
   echo "build ${CUDA_ARCH:-} ${MAKEFLAGS:-}"
   nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
   nvcc --version | tail -n 2
   printf '%s\n' "${stress[@]}" "${small[@]}"
   (cd developers/test/validation && find cases -type f | sort | xargs sha256sum && sha256sum vtest.py)
)
cache=$results/identity-cache/$(echo "$key_text" | sha256sum | cut -c1-16)
[ "${NO_CACHE:-0}" = 1 ] && rm -rf $cache
mkdir -p $cache; echo "$key_text" > $cache/key.txt; touch $cache
ls -dt $results/identity-cache/*/ | tail -n +3 | xargs -r rm -rf

source developers/test/validation/incremental_build.sh

{ echo "reference: $ref ($sha), seeds $seeds, target $target ${MGPU_RANKS:+on $MGPU_RANKS ranks}, cache $cache"; nvidia-smi; nvcc --version; echo "working tree: $(git rev-parse HEAD)"; git status --short; } > $out/system.txt 2>&1

for pr in dp sp; do
   PR=$(echo $pr | tr a-z A-Z)

   cached=0; [ -f $cache/complete-$pr ] && cached=1
   if [ $cached = 0 ]; then
      echo "=== Building reference version ($sha, $pr) ==="
      src=$results/identity-src-$sha
      rm -rf $src; mkdir -p $src
      git archive "$ref" | tar -x -C $src
      if ! make -C $src -j$jobs $target CUDA_PR=-DCUDA_$PR > $out/build-ref-$pr.log 2>&1; then echo "Build failed, see $out/build-ref-$pr.log"; exit 1; fi
      cp $src/$binary $cache/vampire-ref-$pr
      rm -rf $src
   else
      echo "=== Reference version ($sha, $pr) from $cache ==="
   fi
   cp $cache/vampire-ref-$pr ./vampire-cuda-id-ref-$pr

   echo "=== Building working tree version ($pr) ==="
   if ! incremental_build $results/identity-build/$target-$pr $target CUDA_PR=-DCUDA_$PR > $out/build-new-$pr.log 2>&1; then echo "Build failed, see $out/build-new-$pr.log"; exit 1; fi
   cp $results/identity-build/$target-$pr/$binary ./vampire-cuda-id-new-$pr

   for set in "" -stress; do
      rm -rf $results/id-ref-$pr$set
      if [ $cached = 0 ]; then
         echo "=== Validation cases${set:+ with statistics output}, ref version ($pr) ==="
         if [ -z "$set" ]; then $vt run --exe ./vampire-cuda-id-ref-$pr --label id-ref-$pr --precision $pr --seeds $seeds --jobs $gpu_jobs "${ranks[@]}"
         else $vt run --exe ./vampire-cuda-id-ref-$pr --label id-ref-$pr-stress --precision $pr --seeds $seeds --jobs $gpu_jobs "${ranks[@]}" "${extra[@]}"; fi
         rm -rf $cache/id-ref-$pr$set; mv $results/id-ref-$pr$set $cache/
      fi
      ln -s $PWD/$cache/id-ref-$pr$set $results/id-ref-$pr$set
   done
   touch $cache/complete-$pr

   echo "=== Validation cases, new version ($pr) ==="
   $vt run --exe ./vampire-cuda-id-new-$pr --label id-new-$pr --precision $pr --seeds $seeds --jobs $gpu_jobs "${ranks[@]}"
   echo "=== Validation cases with statistics output, new version ($pr) ==="
   $vt run --exe ./vampire-cuda-id-new-$pr --label id-new-$pr-stress --precision $pr --seeds $seeds --jobs $gpu_jobs "${ranks[@]}" "${extra[@]}"

   echo "=== Cases with non-magnetic atoms, with and without gpu:skip-non-magnetic-atoms ($pr) ==="
   $vt run --exe ./vampire-cuda-id-new-$pr --label id-new-$pr-noskip --precision $pr --cases $keep_cases --seeds 1 --jobs $gpu_jobs "${ranks[@]}" "${configs[@]}"
   $vt run --exe ./vampire-cuda-id-new-$pr --label id-new-$pr-skip --precision $pr --cases $keep_cases --seeds 1 --jobs $gpu_jobs "${ranks[@]}" "${configs[@]}" --extra "gpu:skip-non-magnetic-atoms"
done

# skipped atoms change nothing but their own spins
for pr in dp sp; do
   {
      $vt identical --test id-new-$pr-skip --ref id-new-$pr-noskip --cases $keep_cases
      $vt configs --test id-new-$pr-skip --ref id-new-$pr-noskip --cases $keep_cases --magnetic-only
      grep -L "Not integrating" $results/id-new-$pr-skip/*/*/log | sed 's/^/skipping not active: /'
      for d in $results/id-new-$pr-skip/*/*/; do
         run=${d#$results/id-new-$pr-skip/}
         echo "time ${run%/}: $(grep -h 'Simulation run time' $results/id-new-$pr-noskip/$run/stdout | awk '{print $NF}') s without skipping, $(grep -h 'Simulation run time' $d/stdout | awk '{print $NF}') s with"
      done
   } > $out/skip-$pr.txt
   grep -q "FAIL\|not active" $out/skip-$pr.txt && fail_skip=1
done

#------------------------------------------------------------------------------
# byte for byte comparison of every file written by each run
#------------------------------------------------------------------------------
report=$out/identity.txt
fail=${fail_skip:-0}
for pr in dp sp; do
   for set in "" -stress; do
      refdir=$results/id-ref-$pr$set
      for rundir in $(cd $refdir && find . -mindepth 2 -maxdepth 2 -type d | sort); do
         case=$(dirname ${rundir#./})
         # which statistics path the run used (both versions log the same)
         if [ ! -s $refdir/$rundir/output ]; then path="NO OUTPUT"
         elif grep -q "Statistics are calculated on the CPU" $refdir/$rundir/log 2>/dev/null; then path="cpu-stats"
         else path="gpu-stats"; fi
         files=0; diffs=""
         for f in $(cd $refdir/$rundir && find . -type f ! -name stdout ! -name log | sort); do
            files=$((files + 1))
            new=$results/id-new-$pr$set/$rundir/$f
            if [ ! -f $new ] || ! cmp -s <(grep -av '^#' $refdir/$rundir/$f) <(grep -av '^#' $new); then diffs="$diffs ${f#./}"; fi
         done
         # files only the new version wrote
         for f in $(cd $results/id-new-$pr$set/$rundir && find . -type f ! -name stdout ! -name log | sort); do
            [ -f $refdir/$rundir/$f ] || diffs="$diffs ${f#./}(new only)"
         done
         unsupported="GPU acceleration does not support"
         if [ -z "$diffs" ] && [ "$path" != "NO OUTPUT" ]; then status=IDENTICAL
         # features the GPU versions refuse (see run_gpu_suite.sh): both versions must refuse them
         elif [ -z "$diffs" ] && grep -q "$unsupported" $refdir/$rundir/log && grep -q "$unsupported" $results/id-new-$pr$set/$rundir/log; then path="refused"; status="REFUSED by both versions"
         # cases that cannot run on this number of ranks (for example a single spin on several ranks)
         elif [ -z "$diffs" ] && [ -n "${MGPU_RANKS:-}" ] && [ ! -s $results/id-new-$pr$set/$rundir/output ]; then status="FAILED in both versions"
         else status="DIFFERENT:$diffs"; fail=1; fi
         s=${set#-}; printf "%-4s %-7s %-26s %-6s %-10s %3d files  %s\n" $pr "${s:-plain}" $case ${rundir##*/} "$path" $files "$status"
      done
   done
done > $report

#------------------------------------------------------------------------------
# timings: hamr-perf and a 5 nm thermal-perf system (host overhead bound)
#------------------------------------------------------------------------------
for pr in dp sp; do
   for v in ref new; do
      for b in "${benches[@]}"; do
         set=${b% *}; c=${b#* }
         rm -rf $results/id-$v-$pr-$set
         # the reference outputs are cached, so without timings only the working tree version runs
         if [ $v = ref ] && [ "${TIMINGS:-1}" = 0 ] && [ -d $cache/id-ref-$pr-$set ]; then
            ln -s $PWD/$cache/id-ref-$pr-$set $results/id-ref-$pr-$set; continue
         fi
         if [ $set = bench ]; then $vt run --exe ./vampire-cuda-id-$v-$pr --label id-$v-$pr-$set --precision $pr --cases $c > /dev/null
         else $vt run --exe ./vampire-cuda-id-$v-$pr --label id-$v-$pr-$set --precision $pr --cases $c "${small[@]}" > /dev/null; fi
         [ $v = ref ] && { rm -rf $cache/id-ref-$pr-$set; cp -r $results/id-ref-$pr-$set $cache/; }
      done
   done
   for b in "${benches[@]}"; do
      set=${b% *}; c=${b#* }
      for v in ref new; do
         run=run; [ $c = hamr-granular ] && run=seed1
         d=$results/id-$v-$pr-$set/$c/$run
         echo "$pr $c$([ $set = bench-small ] && echo " ($([ -n "${MGPU_RANKS:-}" ] && echo '32 nm, 1 rank' || echo '5 nm'))") $v$([ -L $results/id-$v-$pr-$set ] && echo ' (cached, earlier run)'): $(grep -h 'Simulation run time' $d/stdout | awk '{print $NF}') s, $(grep -c . $d/output) output lines, identical output: $(cmp -s <(grep -av '^#' $results/id-ref-$pr-$set/$c/$run/output) <(grep -av '^#' $d/output) && echo yes || echo NO)"
      done
   done
done | tee $out/benchmarks.txt
grep -q 'identical output: NO' $out/benchmarks.txt && fail=1
rm -f vampire-cuda-id-*

tar czf $results/$name.tar.gz --exclude='*.data' --exclude='*.meta' -C $results $name

echo
echo "=== Summary ==="
for pr in dp sp; do
   for set in plain stress; do
      echo "$pr $set: $(awk -v p=$pr -v s=$set '$1 == p && $2 == s' $report | grep -c IDENTICAL) identical, $(awk -v p=$pr -v s=$set '$1 == p && $2 == s' $report | grep -vc IDENTICAL) different; runs using GPU statistics: $(awk -v p=$pr -v s=$set '$1 == p && $2 == s && $5 == "gpu-stats"' $report | wc -l)"
   done
done
grep -v "IDENTICAL\|REFUSED\|FAILED in both" $report
echo "refused by both versions: $(grep -c REFUSED $report) runs"
[ -n "${MGPU_RANKS:-}" ] && echo "failed in both versions on $MGPU_RANKS ranks: $(grep -c 'FAILED in both' $report) runs ($(grep 'FAILED in both' $report | awk '{print $3}' | sort -u | tr '\n' ' '))"
for pr in dp sp; do
   echo "$pr gpu:skip-non-magnetic-atoms vs without: $(grep -c IDENTICAL $out/skip-$pr.txt) identical, $(grep 'DIFFERENT\|FAIL\|not active' $out/skip-$pr.txt | grep -vc OVERALL) different"
   grep "DIFFERENT\|not active\|^time" $out/skip-$pr.txt | sed 's/^/   /'
done
grep 'identical output: NO' $out/benchmarks.txt
[ $fail = 0 ] && echo "ALL RUNS IDENTICAL" || echo "DIFFERENCES FOUND (see $report)"
echo "Archive: $results/$name.tar.gz"
exit $fail
