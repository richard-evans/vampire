#!/bin/bash
#------------------------------------------------------------------------------
#
#   This file is part of the VAMPIRE open source package under the
#   Free BSD licence (see licence file for details).
#
#------------------------------------------------------------------------------
# Checks that the CPU versions (serial and MPI) and the data converter (vdc)
# built from the working tree behave for an ordinary user exactly like those
# built from a reference commit: every output file must be identical, byte for
# byte. Comment lines starting with '#' (date, process id, path, git hash) in
# text files and the stdout and log files are ignored.
#
#   samples          every sample of samples/ run with both serial versions
#   config           a system with two magnetic and one non-magnetic material
#                    writes spin and cell configuration files in every
#                    combination of config:output-format and config:output-mode,
#                    serial and on RANKS MPI ranks; all files, and everything
#                    each version's vdc makes from them with every vdc option,
#                    must be identical
#   old-files        vdc of the working tree on the files written by the
#                    reference version: its products must equal those of the
#                    reference vdc
#   checkpoints      a 300 K run writes a checkpoint half way (serial and MPI);
#                    continuing it (sim:load-checkpoint=continue) and
#                    restarting from it with every combination of writing and
#                    reading version must give the results of the reference
#                    version reading its own checkpoint
#   window-snapshots (working tree only) a granular HAMR write with
#                    hamr:active-window and spin configuration files: vdc
#                    must give the same products from the window snapshots as
#                    from the full snapshots (hamr:active-window-full-snapshots)
#   warnings         the compiler warnings of each target: the working tree
#                    must not warn in files where the reference does not
#
# Usage (from anywhere):  developers/test/validation/compat_check.sh [ref-commit]
#                         (default: upstream/v8)
#
# Environment variables:
#   JOBS   parallel build jobs (default: number of cores)
#   RANKS  MPI ranks of the parallel runs (default: 2)
#   MPIRUN MPI launcher (default: mpirun --oversubscribe --bind-to none)
#   ONLY   comma separated list of the stages to run (default: all of
#          samples,config,old-files,checkpoints,window-snapshots,warnings)
#   WORK   scratch directory (default: a temporary directory, removed at the
#          end unless KEEP=1); the working tree is always rebuilt, unless
#          REUSE=1 and WORK holds the binaries of an earlier run
#   CACHE  directory for the binaries and build logs of reference commits
#          (default: results/compat-cache), reused by later runs
#
# The summary (PASS and FAIL lines, times) is written to
# developers/test/validation/results/compat-<ref>-<date>/summary.txt
#------------------------------------------------------------------------------
set -u
cd "$(dirname "$0")/../../.."
top=$PWD

ref=${1:-upstream/v8}
jobs=${JOBS:-$(nproc)}
ranks=${RANKS:-2}
only=${ONLY:-samples,config,old-files,checkpoints,window-snapshots,warnings}
mpirun_cmd=${MPIRUN:-mpirun --oversubscribe --bind-to none}
[ "$(id -u)" = 0 ] && mpirun_cmd="$mpirun_cmd --allow-run-as-root"
sha=$(git rev-parse --verify "$ref^{commit}") || exit 1
results=developers/test/validation/results
out=$results/compat-$(echo "$ref" | tr '/' '-')-$(date +%Y%m%d-%H%M)
mkdir -p "$out"; out=$(realpath "$out")
cache=$(realpath -m "${CACHE:-$results/compat-cache}")/$sha
if [ -n "${WORK:-}" ]; then work=$(realpath -m "$WORK"); mkdir -p "$work"; else work=$(mktemp -d); fi
[ "${KEEP:-0}" = 1 ] || [ -n "${WORK:-}" ] || trap 'rm -rf "$work"' EXIT
summary=$out/summary.txt
: > "$summary"
nfail=0

wants(){ case ",$only," in *",$1,"*) return 0;; esac; return 1; }
say(){ echo "$*" | tee -a "$summary"; }
pass(){ say "PASS $*"; }
fail(){ say "FAIL $*"; nfail=$((nfail + 1)); }
note(){ say "note $*"; }
stage_begin(){ stage=$1; t0=$SECONDS; say "=== $1 ==="; }
stage_end(){ say "time $stage: $((SECONDS - t0)) s"; }

#------------------------------------------------------------------------------
# builds: serial, parallel and vdc; binaries and logs kept in a bin directory
#------------------------------------------------------------------------------
build(){ # source-dir bin-dir
   local src=$1 bin=$2
   mkdir -p "$bin"
   ( cd "$src"
     make -j$jobs serial > "$bin/build-serial.log" 2>&1 && cp vampire-serial "$bin/"
     make clean
     make -j$jobs parallel > "$bin/build-parallel.log" 2>&1 && cp vampire-parallel "$bin/"
     make clean
     make vdc > "$bin/build-vdc.log" 2>&1 && cp vdc/vdc "$bin/"
     rm -f vdc/obj/*.o ) > /dev/null 2>&1
   for f in vampire-serial vampire-parallel vdc; do [ -x "$bin/$f" ] || { echo "build of $f failed, see $bin"; exit 1; }; done
}

sources="src hdr obj makefile vdc samples"
t0=$SECONDS
if [ -x "$cache/vdc" ] && [ -x "$cache/vampire-serial" ] && [ -x "$cache/vampire-parallel" ]; then
   echo "reference $ref ($sha) binaries from $cache"
else
   echo "building reference $ref"
   rm -rf "$work/ref-src"; mkdir -p "$work/ref-src"
   git archive "$sha" $sources | tar -x -C "$work/ref-src"
   build "$work/ref-src" "$cache.tmp" && rm -rf "$cache" && mv "$cache.tmp" "$cache"
   rm -rf "$work/ref-src"
fi
if [ "${REUSE:-0}" = 1 ] && [ -x "$work/bin-new/vdc" ] && [ -d "$work/samples" ]; then
   echo "working tree binaries from an earlier run in $work"
else
   echo "building working tree"
   rm -rf "$work/new-src" "$work/bin-new"; mkdir -p "$work/new-src"
   git ls-files -co --exclude-standard -- $sources | tar --ignore-failed-read -cT - | tar -x -C "$work/new-src"
   rm -rf "$work/samples"; cp -r "$work/new-src/samples" "$work/samples"
   build "$work/new-src" "$work/bin-new"
   rm -rf "$work/new-src"
fi
for t in serial parallel vdc; do cp "$cache/build-$t.log" "$out/build-ref-$t.log"; cp "$work/bin-new/build-$t.log" "$out/build-new-$t.log"; done
say "reference: $ref ($(git rev-parse --short "$sha")), working tree: $(git rev-parse --short HEAD)$(git diff --quiet HEAD -- $sources || echo ' + uncommitted changes'), $ranks ranks"
say "time build: $((SECONDS - t0)) s"

declare -A bin=([ref]=$cache [new]=$work/bin-new)

#------------------------------------------------------------------------------
# helpers
#------------------------------------------------------------------------------
# md5 of a file, ignoring comment lines of text files
norm_md5(){ if grep -Iq . "$1" 2>/dev/null; then grep -av '^#' "$1" | md5sum | cut -c1-32; else md5sum < "$1" | cut -c1-32; fi; }

# manifest of the regular files of a directory (symbolic links are skipped), without stdout, log and the given patterns
manifest(){ # dir [excluded name patterns...]
   local d=$1 f p skip; shift
   ( cd "$d" && find . -type f ! -name stdout ! -name log | sort | while read -r f; do
        skip=0; for p in "$@"; do case ${f#./} in $p) skip=1;; esac; done
        [ $skip = 1 ] || echo "$(norm_md5 "$f") ${f#./}"
     done )
}

# names of the files that differ between two manifests (files)
manifest_diff(){ diff "$1" "$2" | grep '^[<>]' | awk '{print $3}' | sort -u | tr '\n' ' '; }

# run a version: dir version kind [ranks]; exit status in dir.status, output in dir/stdout
run_vampire(){
   local dir=$1 v=$2 kind=$3 rc
   if [ "$kind" = serial ]; then ( cd "$dir" && "${bin[$v]}/vampire-serial" > stdout 2>&1 ); rc=$?
   else ( cd "$dir" && $mpirun_cmd -np $ranks "${bin[$v]}/vampire-parallel" > stdout 2>&1 ); rc=$?; fi
   echo $rc > "$dir.status"
   return $rc
}

# compare the files of two run directories; prints PASS/FAIL with the label
compare_runs(){ # label dir-ref dir-new [excluded patterns...]
   local label=$1 a=$2 b=$3; shift 3
   manifest "$a" "$@" > "$a.man"; manifest "$b" "$@" > "$b.man"
   local d; d=$(manifest_diff "$a.man" "$b.man")
   local n; n=$(wc -l < "$a.man")
   if [ "$n" = 0 ]; then fail "$label: no files written"
   elif [ -z "$d" ]; then pass "$label: $n files identical"; else fail "$label: differing files: $d"; fi
}

#------------------------------------------------------------------------------
# 1. samples
#------------------------------------------------------------------------------
if wants samples; then
   stage_begin samples
   for input in $(cd "$work/samples" && find . -name input | sort); do
      s=$(dirname "${input#./}")
      for v in ref new; do
         d=$work/samples-run/$v/$s; mkdir -p "$d"; cp -r "$work/samples/$s/." "$d/"
         if [ -f "$d/unit-cell-creator.py" ]; then ( cd "$d" && python3 unit-cell-creator.py > /dev/null 2>&1 ); fi
         steps=$(sed -n 's/^ *sim:total-time-steps *= *\([0-9]*\).*/\1/p' "$d/input" | tail -n 1)
         if [ -n "$steps" ] && [ "$steps" -gt 50000 ]; then
            sed -i "s/^\( *sim:total-time-steps *=\).*/\1 50000/" "$d/input"
            [ $v = ref ] && note "sample $s: total time steps shortened from $steps to 50000 for both versions"
         fi
         run_vampire "$d" $v serial
      done
      ra=$(cat "$work/samples-run/ref/$s.status"); rb=$(cat "$work/samples-run/new/$s.status")
      if [ "$ra" != 0 ] || [ "$rb" != 0 ]; then
         if [ "$ra" = "$rb" ]; then note "sample $s: both versions exit with status $ra, no comparison"; else fail "sample $s: exit status $ra (reference) and $rb (working tree)"; fi
         continue
      fi
      compare_runs "sample $s" "$work/samples-run/ref/$s" "$work/samples-run/new/$s"
   done
   rm -rf "$work/samples-run"
   stage_end
fi

#------------------------------------------------------------------------------
# small test system shared by the config and checkpoints stages
#------------------------------------------------------------------------------
common=$work/common; mkdir -p "$common"
cat > "$common/m.mat" << 'EOF'
material:num-materials=3
material[1]:material-name=A
material[1]:damping-constant=0.1
material[1]:exchange-matrix[1]=11.2e-21
material[1]:exchange-matrix[2]=-2.0e-21
material[1]:exchange-matrix[3]=0
material[1]:atomic-spin-moment=1.72 !muB
material[1]:second-order-uniaxial-anisotropy-constant=1.0e-23
material[1]:minimum-height=0.0
material[1]:maximum-height=0.4
material[1]:initial-spin-direction=random
material[2]:material-name=B
material[2]:damping-constant=0.3
material[2]:exchange-matrix[1]=-2.0e-21
material[2]:exchange-matrix[2]=6.0e-21
material[2]:exchange-matrix[3]=0
material[2]:atomic-spin-moment=3.0 !muB
material[2]:second-order-uniaxial-anisotropy-constant=5.0e-23
material[2]:minimum-height=0.4
material[2]:maximum-height=0.8
material[2]:initial-spin-direction=random
material[3]:material-name=Cap
material[3]:non-magnetic=remove
material[3]:minimum-height=0.8
material[3]:maximum-height=1.0
material[3]:exchange-matrix[1]=0
material[3]:exchange-matrix[2]=0
material[3]:exchange-matrix[3]=0
material[3]:atomic-spin-moment=1.0 !muB
EOF
cat > "$common/base" << 'EOF'
create:crystal-structure=sc
dimensions:unit-cell-size = 3.54 !A
dimensions:system-size-x = 6 !nm
dimensions:system-size-y = 6 !nm
dimensions:system-size-z = 4 !nm
material:file=m.mat
sim:temperature=0.0
sim:time-steps-increment=1
sim:total-time-steps=300
sim:time-step=1.0E-15
sim:applied-field-strength=0.5 !T
sim:applied-field-unit-vector=1,0,1
sim:program=time-series
sim:integrator=llg-heun
output:real-time
output:magnetisation
output:material-magnetisation
output:precision=16
config:atoms
config:atoms-output-rate=100
config:macro-cells
config:macro-cells-output-rate=100
EOF

#------------------------------------------------------------------------------
# 2. config and 3. old-files: configuration output, vdc
#------------------------------------------------------------------------------
vdc_options=("" "--xyz" "--vtk" "--povray" "--text" "--cells" "--ssc" "--grains" "--povray-cells" "--povray-sticks"
             "--cell-size 10 10 10" "--cells --cell-size 10 10 10" "--povray --povray-sticks" "--xyz --grains")

# all products of vdc (version v) on the files of a run directory: manifests in dir.vdc-<v>/<n>, exit codes in .rc
vdc_all(){ # version run-dir out-dir
   local v=$1 src=$2 dst=$3 n=0 opt p
   mkdir -p "$dst"
   for opt in "${vdc_options[@]}"; do
      n=$((n + 1)); p=$dst/pd$n; rm -rf "$p"; mkdir "$p"
      for f in "$src"/*; do case $(basename "$f") in stdout|log) ;; *) [ -f "$f" ] && ln -s "$(realpath "$f")" "$p/";; esac; done
      ( cd "$p" && "${bin[$v]}/vdc" $opt > "$dst/out$n" 2>&1 ); echo $? > "$dst/rc$n"
      manifest "$p" > "$dst/man$n"
      # an error message is part of the behaviour; progress output is not
      if [ "$(cat "$dst/rc$n")" != 0 ]; then cp "$dst/out$n" "$dst/err$n"; else : > "$dst/err$n"; fi
      rm -rf "$p"
   done
}

vdc_compare(){ # label dir-a dir-b
   local label=$1 a=$2 b=$3 n=0 opt bad=""
   for opt in "${vdc_options[@]}"; do
      n=$((n + 1))
      if [ "$(cat "$a/rc$n")" != "$(cat "$b/rc$n")" ] || ! cmp -s "$a/man$n" "$b/man$n" || ! cmp -s "$a/err$n" "$b/err$n"; then
         bad="$bad [vdc $opt: exit $(cat "$a/rc$n") / $(cat "$b/rc$n"); files $(manifest_diff "$a/man$n" "$b/man$n")]"
      fi
   done
   if [ -z "$bad" ]; then pass "$label: vdc products identical for ${#vdc_options[@]} option sets ($(cat "$a"/man* | wc -l) product files)"; else fail "$label:$bad"; fi
}

if wants config || wants old-files; then
   stage_begin config
   cfg=$work/config; rm -rf "$cfg"; mkdir -p "$cfg"
   # combinations: name format mode nodes
   modes=()
   for fmt in text binary; do
      for mode in legacy mpi-io file-per-process; do modes+=("$fmt-$mode $fmt $mode 1"); done
      modes+=("$fmt-file-per-node-1 $fmt file-per-node 1")
      [ $ranks -ge 2 ] && modes+=("$fmt-file-per-node-2 $fmt file-per-node 2")
   done
   for kind in serial parallel; do
      for m in "${modes[@]}"; do
         set -- $m; name=$1
         [ $kind = serial ] && [ "$4" = 2 ] && continue
         for v in ref new; do
            d=$cfg/$v/$kind/$name; mkdir -p "$d"; cp "$common/m.mat" "$d/"
            { cat "$common/base"; echo "config:output-format=$2"; echo "config:output-mode=$3"; echo "config:output-nodes=$4"; } > "$d/input"
            run_vampire "$d" $v $kind
         done
         a=$cfg/ref/$kind/$name; b=$cfg/new/$kind/$name
         ra=$(cat "$a.status"); rb=$(cat "$b.status")
         if [ "$ra" != 0 ] || [ "$rb" != 0 ]; then
            if [ "$ra" = "$rb" ]; then
               note "$kind $name: both versions refuse the combination ($(grep -ai 'error' "$a/stdout" | head -n 1))"
            else fail "$kind $name: exit status $ra (reference) and $rb (working tree)"; fi
            touch "$a.refused"; continue
         fi
         compare_runs "config $kind $name" "$a" "$b"
         if wants config; then
            if [ -f "$a/atoms-coords.meta" ]; then
               vdc_all ref "$a" "$a.vdc-ref"; vdc_all new "$b" "$b.vdc-new"
               vdc_compare "config $kind $name" "$a.vdc-ref" "$b.vdc-new"
               rm -rf "$b.vdc-new"
            else
               # legacy files are not vdc input: both versions must stop with the same error
               rc_a=$( cd "$a" && "${bin[ref]}/vdc" --xyz > "$a.vdc-out" 2>&1; echo $? ); rc_b=$( cd "$b" && "${bin[new]}/vdc" --xyz > "$b.vdc-out" 2>&1; echo $? )
               rm -f "$a/crystal.xyz" "$b/crystal.xyz"
               if [ "$rc_a" = "$rc_b" ] && [ "$rc_a" != 0 ] && cmp -s "$a.vdc-out" "$b.vdc-out"; then pass "config $kind $name: vdc stops with the same error in both versions"
               else fail "config $kind $name: vdc on legacy files: exit $rc_a / $rc_b, messages $(cmp -s "$a.vdc-out" "$b.vdc-out" && echo same || echo different)"; fi
            fi
         fi
      done
   done
   stage_end

   if wants old-files; then
      stage_begin old-files
      for kind in serial parallel; do
         for m in "${modes[@]}"; do
            set -- $m; name=$1; a=$cfg/ref/$kind/$name
            [ -f "$a.refused" ] || [ ! -f "$a/atoms-coords.meta" ] && continue
            [ -d "$a.vdc-ref" ] || vdc_all ref "$a" "$a.vdc-ref"
            vdc_all new "$a" "$a.vdc-new"
            vdc_compare "old files $kind $name" "$a.vdc-ref" "$a.vdc-new"
            rm -rf "$a.vdc-new"
         done
      done
      stage_end
   fi
   rm -rf "$cfg"
fi

#------------------------------------------------------------------------------
# 4. checkpoints
#------------------------------------------------------------------------------
if wants checkpoints; then
   stage_begin checkpoints
   ck=$work/checkpoints; rm -rf "$ck"; mkdir -p "$ck"
   cp "$common/m.mat" "$ck/m.mat"
   grep -v '^config\|^sim:total-time\|^sim:temperature' "$common/base" > "$ck/base"
   ckopts=("sim:temperature=300" "sim:integrator-random-seed=4242" "output:output-rate=10" "config:atoms=end")
   for kind in serial parallel; do
      # writers: 150 of 300 steps, checkpoint at the end
      for w in ref new; do
         d=$ck/w-$w-$kind; mkdir -p "$d"; cp "$ck/m.mat" "$d/"
         { cat "$ck/base"; printf '%s\n' "${ckopts[@]}" "sim:total-time-steps=150" "sim:save-checkpoint=end"; } > "$d/input"
         run_vampire "$d" $w $kind || fail "checkpoint writer $w $kind: exit status $(cat "$d.status")"
      done
      compare_runs "checkpoint writers $kind (files other than the checkpoint)" "$ck/w-ref-$kind" "$ck/w-new-$kind" 'vampire*.chk'
      # a checkpoint of the working tree may hold more data, but must begin with the whole reference checkpoint
      for c in "$ck/w-ref-$kind"/vampire*.chk; do
         c=$(basename "$c"); n=$(stat -c %s "$ck/w-ref-$kind/$c"); m=$(stat -c %s "$ck/w-new-$kind/$c")
         if [ "$n" = "$m" ] && cmp -s "$ck/w-ref-$kind/$c" "$ck/w-new-$kind/$c"; then pass "checkpoint $kind $c: identical"
         elif [ "$m" -gt "$n" ] && cmp -s -n "$n" "$ck/w-ref-$kind/$c" "$ck/w-new-$kind/$c"; then note "checkpoint $kind $c: the working tree file is the reference file with $((m - n)) bytes appended"
         else fail "checkpoint $kind $c: differs from the reference file ($n and $m bytes)"; fi
      done
      # readers: continue to 300 steps, and restart (300 steps from the saved spins)
      for mode in continue restart; do
         for w in ref new; do
            for r in ref new; do
               d=$ck/r-$mode-$w-$r-$kind; mkdir -p "$d"; cp "$ck/m.mat" "$ck"/w-$w-$kind/vampire*.chk "$d/"
               { cat "$ck/base"; printf '%s\n' "${ckopts[@]}" "sim:total-time-steps=300" "sim:load-checkpoint=$mode"; } > "$d/input"
               run_vampire "$d" $r $kind || fail "checkpoint $mode $kind, written by $w, read by $r: exit status $(cat "$d.status")"
            done
         done
         base=$ck/r-$mode-ref-ref-$kind
         for w in ref new; do
            for r in ref new; do
               [ $w = ref ] && [ $r = ref ] && continue
               compare_runs "checkpoint $mode $kind, written by $w, read by $r, against reference -> reference" "$base" "$ck/r-$mode-$w-$r-$kind" 'vampire*.chk'
            done
         done
      done
      # the output file of a continued run starts at the checkpoint time, so it must differ from the restarted run
      if cmp -s <(grep -av '^#' "$ck/r-continue-ref-ref-$kind/output") <(grep -av '^#' "$ck/r-restart-ref-ref-$kind/output"); then
         fail "checkpoint $kind: continue and restart give the same output, the check does not distinguish them"
      fi
   done
   rm -rf "$ck"
   stage_end
fi

#------------------------------------------------------------------------------
# 5. window snapshots (working tree only)
#------------------------------------------------------------------------------
if wants window-snapshots; then
   stage_begin window-snapshots
   ws=$work/window; rm -rf "$ws"; mkdir -p "$ws"
   hcase=$top/developers/test/validation/cases/hamr-granular
   sed -e 's/^material\[4\]:density = 1.0/material[4]:density = 0.01/' \
       -e 's/initial-spin-direction = 0.3,0,1/initial-spin-direction = 0,0,1/' \
       -e 's/^material\[5\]:minimum-height=0.5/material[5]:minimum-height=0.85/' \
       "$hcase/FePt_sensor.mat" > "$ws/FePt_sensor.mat"
   sed -e 's/system-size-x = 8 !nm/system-size-x = 20 !nm/' \
       -e 's/system-size-z = 20 !nm/system-size-z = 4 !nm/' \
       -e 's/grain-height = 10 !nm/grain-height = 2 !nm/' \
       -e 's/laser-FWHM-\([xy]\) = 5/laser-FWHM-\1 = 8/' \
       -e 's/head-speed = 5000/head-speed = 3000/' \
       -e 's/bit-size = 4 !nm/bit-size = 6 !nm/' \
       -e 's/number-of-bits = 1$/number-of-bits = 2/' \
       -e 's/equilibration-time-steps = 0/equilibration-time-steps = 1000/' \
       -e 's/output-rate = 500/output-rate = 1000/' \
       "$hcase/input" > "$ws/input.base"
   for kind in serial parallel; do
      for v in window full; do
         d=$ws/$v-$kind; mkdir -p "$d"; cp "$ws/FePt_sensor.mat" "$d/"
         { cat "$ws/input.base"; echo "sim:integrator-random-seed=1"; echo "hamr:active-window"
           [ $v = full ] && echo "hamr:active-window-full-snapshots"
           echo "config:atoms"; echo "config:output-format=binary"; echo "config:atoms-output-rate=5000"; } > "$d/input"
         run_vampire "$d" new $kind || fail "window snapshots $kind $v: exit status $(cat "$d.status")"
      done
      nwin=$(grep -l 'Spin data: window' "$ws/window-$kind"/spins-*.meta 2> /dev/null | wc -l)
      nfull=$(grep -l 'Spin data: window' "$ws/full-$kind"/spins-*.meta 2> /dev/null | wc -l)
      nsnap=$(ls "$ws/window-$kind"/spins-*.meta 2> /dev/null | wc -l)
      if [ "$nwin" -gt 0 ] && [ "$nfull" = 0 ] && [ "$nsnap" -ge 2 ]; then pass "window snapshots $kind: $nwin of $nsnap snapshots are window snapshots, none with full snapshots"
      else fail "window snapshots $kind: $nwin window snapshots of $nsnap (with full snapshots: $nfull)"; fi
      compare_runs "window snapshots $kind output (window against full snapshots)" "$ws/window-$kind" "$ws/full-$kind" 'spins-*' 'input'
      vdc_all new "$ws/window-$kind" "$ws/window-$kind.vdc"; vdc_all new "$ws/full-$kind" "$ws/full-$kind.vdc"
      n=0; bad=""
      for opt in "${vdc_options[@]}"; do
         n=$((n + 1))
         case "$opt" in "--xyz"|"--vtk"|"--text"|"--povray") ;; *) continue;; esac
         if [ "$(cat "$ws/window-$kind.vdc/rc$n")" != 0 ] || ! cmp -s "$ws/window-$kind.vdc/man$n" "$ws/full-$kind.vdc/man$n"; then
            bad="$bad [vdc $opt: exit $(cat "$ws/window-$kind.vdc/rc$n") / $(cat "$ws/full-$kind.vdc/rc$n"); files $(manifest_diff "$ws/window-$kind.vdc/man$n" "$ws/full-$kind.vdc/man$n")]"
         fi
      done
      if [ -z "$bad" ]; then pass "window snapshots $kind: vdc --xyz --vtk --text --povray identical from window and full snapshots"; else fail "window snapshots $kind:$bad"; fi
   done
   rm -rf "$ws"
   stage_end
fi

#------------------------------------------------------------------------------
# 6. warnings
#------------------------------------------------------------------------------
if wants warnings; then
   stage_begin warnings
   for t in serial parallel vdc; do
      files(){ grep -ao '^[^ :]*:[0-9]*:[0-9]*: warning' "$1" | cut -d: -f1 | sort | uniq -c | awk '{print $2, $1}'; }
      files "$cache/build-$t.log" > "$work/w-ref"; files "$out/build-new-$t.log" > "$work/w-new"
      nr=$(awk '{s += $2} END {print s + 0}' "$work/w-ref"); nn=$(awk '{s += $2} END {print s + 0}' "$work/w-new")
      newonly=$(awk 'NR == FNR {r[$1]; next} !($1 in r) {printf "%s(%d) ", $1, $2}' "$work/w-ref" "$work/w-new")
      if [ -z "$newonly" ]; then pass "warnings $t: $nr (reference), $nn (working tree), none in files the reference does not warn about"
      else fail "warnings $t: $nr (reference), $nn (working tree), new warnings in: $newonly"; fi
   done
   stage_end
fi

say "total time: $SECONDS s, failures: $nfail"
echo "summary: $summary"
[ $nfail = 0 ]
