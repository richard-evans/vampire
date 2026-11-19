#!/bin/bash
#------------------------------------------------------------------------------
#
#   This file is part of the VAMPIRE open source package under the
#   Free BSD licence (see licence file for details).
#
#------------------------------------------------------------------------------
# Compares the CUDA version of vampire built from a reference commit with the
# CUDA version built from the working tree, for changes that are expected to
# give statistically equivalent rather than identical results.
#
#   1. builds both versions in the requested precision
#   2. runs every statistical and analytic-thermal case with both versions,
#      using more seeds than the validation suite for a sharper comparison
#   3. compares new vs reference version and each against the CPU references,
#      with a pooled z-score summary (vtest.py compare --summary)
#   4. times the hamr-perf and thermal-perf benchmarks with both versions
#
# Usage (from anywhere):  developers/test/validation/ab_compare.sh [ref-commit] [dp|sp] [seeds]
#                         (defaults: HEAD, dp, 32)
#
# Environment variables:
#   CASES   comma separated cases to compare (default: all finite temperature LLG cases)
#   JOBS    parallel build jobs (default: number of cores)
#
# Results are packed into developers/test/validation/results/ab-<ref>-<pr>-<host>-<date>.tar.gz
#
# Note: the reference version is built from an exported copy of the commit;
# the working tree build runs 'make clean' before and after.
#------------------------------------------------------------------------------
set -u
cd "$(dirname "$0")/../../.."

ref=${1:-HEAD}
pr=${2:-dp}
seeds=${3:-32}
jobs=${JOBS:-$(nproc)}
sha=$(git rev-parse --short "$ref") || exit 1
PR=$(echo $pr | tr a-z A-Z)

name=ab-$sha-$pr-$(hostname -s)-$(date +%Y%m%d-%H%M)
results=developers/test/validation/results
out=$results/$name
mkdir -p $out
vt="python3 developers/test/validation/vtest.py"
# finite temperature LLG cases (Monte Carlo cases use a separate GPU code path)
statcases=$(python3 -c "
import json, os
d = 'developers/test/validation/cases'
def llg_thermal(c):
   case = json.load(open(os.path.join(d, c, 'case.json')))
   return case['kind'] in ('statistical', 'analytic') and case.get('model') != 'precession' \
      and 'monte-carlo' not in open(os.path.join(d, c, 'input')).read()
print(','.join(sorted(c for c in os.listdir(d) if llg_thermal(c))))")
statcases=${CASES:-$statcases}

{ echo "reference: $ref ($sha), precision $pr, seeds $seeds"; nvidia-smi; nvcc --version; echo "working tree: $(git rev-parse HEAD)"; git status --short; } > $out/system.txt 2>&1

echo "=== Building reference version ($sha, $pr) ==="
src=$results/ab-src-$sha
rm -rf $src; mkdir -p $src
git archive "$ref" | tar -x -C $src
if ! make -C $src -j$jobs cuda CUDA_PR=-DCUDA_$PR > $out/build-ref.log 2>&1; then echo "Build failed, see $out/build-ref.log"; exit 1; fi
cp $src/vampire-cuda ./vampire-cuda-ab-ref
rm -rf $src

echo "=== Building working tree version ($pr) ==="
make clean > /dev/null
if ! make -j$jobs cuda CUDA_PR=-DCUDA_$PR > $out/build-new.log 2>&1; then echo "Build failed, see $out/build-new.log"; exit 1; fi
mv vampire-cuda ./vampire-cuda-ab-new
make clean > /dev/null

for v in ref new; do
   echo "=== Statistical cases, $v version, $seeds seeds ==="
   $vt run --exe ./vampire-cuda-ab-$v --label ab-$v-$pr --precision $pr --cases $statcases --seeds $seeds
done

$vt compare --test ab-new-$pr --ref ab-ref-$pr --cases $statcases --verbose --summary > $out/compare-new-vs-ref.txt
$vt compare --test ab-new-$pr --ref cpu --cases $statcases --verbose --summary > $out/compare-new-vs-cpu.txt
$vt compare --test ab-ref-$pr --ref cpu --cases $statcases --verbose --summary > $out/compare-ref-vs-cpu.txt

for v in ref new; do
   echo "=== Benchmarks, $v version ==="
   $vt run --exe ./vampire-cuda-ab-$v --label ab-$v-$pr-bench --precision $pr --cases hamr-perf,thermal-perf
done
for b in hamr-perf thermal-perf; do
   for v in ref new; do
      echo "$b $v: $(grep -h 'Simulation run time' $results/ab-$v-$pr-bench/$b/run/stdout | awk '{print $NF}') s"
   done
done | tee $out/benchmarks.txt
rm -f vampire-cuda-ab-ref vampire-cuda-ab-new

tar czf $results/$name.tar.gz --exclude='*.data' --exclude='*.meta' -C $results $name ab-ref-$pr ab-new-$pr ab-ref-$pr-bench ab-new-$pr-bench

echo
echo "=== Summary ==="
for f in $out/compare-*.txt; do
   echo "$(basename $f):"; grep -E "^   all|fraction|OVERALL" $f
done
cat $out/benchmarks.txt
echo "Archive: $results/$name.tar.gz"
