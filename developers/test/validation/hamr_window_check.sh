#!/bin/bash
#------------------------------------------------------------------------------
#
#   This file is part of the VAMPIRE open source package under the
#   Free BSD licence (see licence file for details).
#
#------------------------------------------------------------------------------
# Checks hamr:active-window on HAMR writes of 9 bits on a granular FePt film
# (2 nm grain height, laser FWHM 8 nm, 1000 m/s): one track along a 60 x 8 nm
# film, or three tracks of 3 bits on a 20 x 80 nm film (wide film); with FILM=small
# every stage uses one track of 4 bits on a 30 x 8 nm film (minutes on one CPU):
#
#   1. identity: on the wide film of magnetic grains only, the window with
#      hamr:active-window-temperature = 0 (every grain integrated) gives output
#      identical to the run without the window
#   2. statistics: on the 60 x 8 nm film (the wide film with FILM=wide) with a
#      non-magnetic fill and sensor layer, SEEDS
#      seeds with and without the window (threshold 1 K): z-scores of the
#      seed-mean system magnetisation at every output and of the final
#      magnetisation of every grain; FAIL if any |z| exceeds ZMAX (default 4)
#      and the statistics of one window run calculated from the spins on the
#      host agree with the device statistics (largest difference below TOL,
#      default 1e-10)
#   3. window snapshots: one window run with spin snapshots every 5 ps, written
#      as window snapshots and in full: the output is identical, and every
#      snapshot rebuilt from the window snapshots is identical to the full one;
#      vdc (VDC, default vdc/vdc next to the binary, skipped if there is none)
#      gives identical files from both
#   4. the run times with and without the window
#   5. with PERF=<seconds>: the speed on a 300 x 20 x 20 nm film (4.3 million
#      atoms), each run stopped after PERF seconds (600 is enough on a small GPU)
#
#    developers/test/validation/hamr_window_check.sh [vampire binary]
#
# The binary defaults to ./vampire-serial (CUDA binaries run on the GPU).
# SEEDS defaults to 16, JOBS (runs at a time) to 1. With RANKS=N every run uses
# mpirun -np N, or $MPIRUN -np N (e.g. MPIRUN="mpirun --bind-to none" for
# several runs at a time) (MPI and multi-GPU versions; the window then splits the film
# across the track, while runs without it use the usual decomposition). The
# multi-GPU version runs with gpu:decomposition-independent-order, so its
# results do not depend on the decomposition. The 20 x 80 nm film splits across
# the track at 2, 4 and 6 processes without the window too, so the identity
# holds for the MPI version. Results go to a temporary directory, which is
# removed unless KEEP=1, or to the (new, empty) directory WORK, which is kept.
#------------------------------------------------------------------------------
set -eu

exe=$(realpath "${1:-./vampire-serial}")
seeds=${SEEDS:-16}
jobs=${JOBS:-1}
zmax=${ZMAX:-4}
film=${FILM:-narrow}
# films of the identity and snapshot stages (all small with FILM=small)
id_film=wide; snap_film=narrow
[ "$film" = small ] && id_film=small && snap_film=small
mpi=""
[ -n "${RANKS:-}" ] && mpi="${MPIRUN:-mpirun} -np $RANKS"
mgpu=0
case $(basename "$exe") in *multigpu*) mgpu=1;; esac
here=$(dirname "$(realpath "$0")")
if [ -n "${WORK:-}" ]; then
   work=$(realpath -m "$WORK"); mkdir -p "$work"; KEEP=1
else
   work=$(mktemp -d)
fi
[ "${KEEP:-0}" = 1 ] && echo "results in $work" || trap 'rm -rf "$work"' EXIT

# material file of the granular HAMR case: low density fill, thin sensor layer,
# and a version with the three FePt materials only
sed -e 's/^material\[4\]:density = 1.0/material[4]:density = 0.01/' \
    -e 's/initial-spin-direction = 0.3,0,1/initial-spin-direction = 0,0,1/' \
    -e 's/^material\[5\]:minimum-height=0.5/material[5]:minimum-height=0.85/' \
    "$here/cases/hamr-granular/FePt_sensor.mat" > "$work/full.mat"
grep -v '^material\[[45]\]' "$work/full.mat" | sed 's/^material:num-materials *= *5/material:num-materials=3/' > "$work/grains.mat"
grep -q '^material:num-materials=3' "$work/grains.mat"

run(){ # name material film (narrow or wide) seed [extra input lines]
   local dir="$work/$1" mat=$2 size=60x8 bits=9 seed=$4
   [ "$3" = wide ] && size=20x80
   [ "$3" = small ] && size=30x8 && bits=4
   shift 4
   mkdir "$dir"
   cp "$work/$mat" "$dir/FePt_sensor.mat"
   sed -e "s/system-size-x = 8 !nm/system-size-x = ${size%x*} !nm/" \
       -e "s/system-size-y = 8 !nm/system-size-y = ${size#*x} !nm/" \
       -e 's/system-size-z = 20 !nm/system-size-z = 4 !nm/' \
       -e 's/grain-height = 10 !nm/grain-height = 2 !nm/' \
       -e 's/laser-FWHM-\([xy]\) = 5/laser-FWHM-\1 = 8/' \
       -e 's/head-speed = 5000/head-speed = 1000/' \
       -e 's/bit-size = 4 !nm/bit-size = 6 !nm/' \
       -e "s/number-of-bits = 1\$/number-of-bits = $bits/" \
       -e 's/equilibration-time-steps = 0/equilibration-time-steps = 20000/' \
       -e 's/output-rate = 500/output-rate = 1000/' \
       "$here/cases/hamr-granular/input" > "$dir/input"
   echo "sim:integrator-random-seed = $seed" >> "$dir/input"
   [ "$mgpu" = 1 ] && echo "gpu:decomposition-independent-order" >> "$dir/input"
   for line in "$@"; do echo "$line" >> "$dir/input"; done
   local start=$(date +%s.%N)
   (cd "$dir" && $mpi "$exe" > stdout 2>&1) || { tail -n 20 "$dir/stdout"; return 1; }
   echo "$(date +%s.%N) - $start" | bc > "$dir/seconds"
}
export -f run
export work exe here mpi mgpu

fail=0

# 1. identity
run id-off grains.mat $id_film 1
run id-window grains.mat $id_film 1 "hamr:active-window" "hamr:active-window-temperature = 0"
for f in output grain.txt; do
   if cmp -s <(grep -v '^#' "$work/id-off/$f") <(grep -v '^#' "$work/id-window/$f"); then
      echo "PASS identity: $f identical with every grain in the window"
   else
      echo "FAIL identity: $f differs with every grain in the window"; fail=1
   fi
done

# 2. statistics
for s in $(seq 1 "$seeds"); do
   echo "run off-$s full.mat $film $s"
   echo "run window-$s full.mat $film $s hamr:active-window"
done | xargs -P "$jobs" -I{} bash -c '{}'

python3 - "$work" "$seeds" "$zmax" << 'EOF' || fail=1
import sys, numpy as np
work, seeds, zmax = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
rows = lambda f: np.array([l.split() for l in open(f) if l[:1].isdigit()], float)
def z(a, b):
    se = np.sqrt(a.var(0, ddof=1)/len(a) + b.var(0, ddof=1)/len(b))
    return np.where(se > 0, (b.mean(0) - a.mean(0))/np.where(se > 0, se, 1), 0)
o = {m: np.array([rows(f'{work}/{m}-{s}/output') for s in range(1, seeds+1)]) for m in ('off', 'window')}
mz = {m: o[m][:, :, 5]*o[m][:, :, 6] for m in o}
zs = z(mz['off'], mz['window'])
# grain file: time, temperature, then mx, my, mz, |m| of every grain
g = {m: np.array([rows(f'{work}/{m}-{s}/grain.txt')[-1] for s in range(1, seeds+1)]) for m in o}
ngrains = int(open(f'{work}/off-1/grain-coordinates.txt').read().split('Number of grains:')[1].split()[0])
gmz = {m: np.array([g[m][:, 2+4*i+2]*g[m][:, 2+4*i+3] for i in range(ngrains)]).T for m in g}
zg = z(gmz['off'], gmz['window'])
print(f'system mz at {len(zs)} outputs: max |z| {np.abs(zs).max():.2f}; final {mz["off"][:, -1].mean():+.4f} (off) {mz["window"][:, -1].mean():+.4f} (window)')
w = np.abs(zs).argmax()
print(f'   largest at t = {o["off"][0, w, 0]:.3g} s: off {mz["off"][:, w].mean():+.4f} sd {mz["off"][:, w].std(ddof=1):.4f}, window {mz["window"][:, w].mean():+.4f} sd {mz["window"][:, w].std(ddof=1):.4f}')
print(f'final mz of {ngrains} grains: max |z| {np.abs(zg).max():.2f}, sum z^2 {np.sum(zg**2):.1f}')
worst = max(np.abs(zs).max(), np.abs(zg).max())
if seeds < 8:
    # with so few seeds the standard deviations, and so the z-scores, are too uncertain to judge
    print(f'NOTE statistics: largest |z| {worst:.2f} with {seeds} seeds, not judged (at least 8 needed)')
    sys.exit(0)
print(('PASS' if worst <= zmax else 'FAIL') + f' statistics: largest |z| {worst:.2f} (limit {zmax}, {seeds} seeds)')
sys.exit(int(worst > zmax))
EOF

# statistics of a window run calculated from the spins on the host: the same
# trajectory, so the output must agree with the device statistics to rounding
run window-cpustats-1 full.mat $film 1 "hamr:active-window" "gpu:calculate-statistics-on-cpu=true"
python3 - "$work" "${TOL:-1e-10}" << 'EOF' || fail=1
import sys, numpy as np
work, tol = sys.argv[1], float(sys.argv[2])
rows = lambda f: np.array([l.split() for l in open(f) if l[:1].isdigit()], float)
d = max(np.abs(rows(f'{work}/window-1/{f}') - rows(f'{work}/window-cpustats-1/{f}')).max() for f in ('output', 'grain.txt'))
print(('PASS' if d <= tol else 'FAIL') + f' window statistics: device and host statistics differ by at most {d:.2g} (limit {tol})')
sys.exit(int(d > tol))
EOF

# 3. window snapshots
snapshots=(config:atoms "config:output-format = binary" "config:atoms-output-rate = 5000")
run snap-window-1 full.mat $snap_film 1 "hamr:active-window" "${snapshots[@]}"
run snap-full-1 full.mat $snap_film 1 "hamr:active-window" "hamr:active-window-full-snapshots" "${snapshots[@]}"
python3 - "$work" << 'EOF' || fail=1
import sys, os, re, numpy as np
work = sys.argv[1]
rows = lambda f: [l for l in open(f) if not l.startswith('#')]
def series(run):
    state, i = None, 0
    while os.path.exists(meta := f'{run}/spins-{i:08d}.meta'):
        text = open(meta).read()
        data = f'{run}/' + re.search(r'^(spins-\S+\.data)$', text, re.M).group(1)
        with open(data, 'rb') as f:
            if 'Spin data: window' in text:
                n, r = np.fromfile(f, np.uint64, 2)
                ranges = np.fromfile(f, np.uint64, 2 * int(r)).astype(int).reshape(-1, 2)
                s = np.fromfile(f, np.float64).reshape(-1, 3)
                pos = np.concatenate([np.arange(a, a + c) for a, c in ranges]) if len(ranges) else []
                state[pos] = s
            else:
                n = int(np.fromfile(f, np.uint64, 1)[0])
                state = np.fromfile(f, np.float64).reshape(n, 3)
        yield 'Spin data: window' in text, os.path.getsize(data), state.copy()
        i += 1
ok = all(rows(f'{work}/snap-window-1/{f}') == rows(f'{work}/snap-full-1/{f}') for f in ('output', 'grain.txt'))
w, f = list(series(f'{work}/snap-window-1')), list(series(f'{work}/snap-full-1'))
same = len(w) == len(f) and all(np.array_equal(a[2], b[2]) for a, b in zip(w, f))
nw = sum(a[0] for a in w)
print(('PASS' if ok and same and nw > 0 else 'FAIL') + f' window snapshots: {len(w)} snapshots ({nw} window), output '
      f'{"identical" if ok else "differs"}, rebuilt snapshots {"identical" if same else "differ"}; spin data '
      f'{sum(a[1] for a in w) / 1e6:.1f} MB instead of {sum(b[1] for b in f) / 1e6:.1f} MB')
sys.exit(int(not (ok and same and nw > 0)))
EOF

# vdc on the window and the full snapshots: every file it writes must be identical
vdc=$(realpath -m "${VDC:-$(dirname "$exe")/vdc/vdc}")
if [ -f "$vdc" ] && [ -x "$vdc" ]; then
   for d in snap-window-1 snap-full-1; do
      ls "$work/$d" > "$work/$d.before"
      (cd "$work/$d" && "$vdc" --xyz --vtk --text > vdc.log 2>&1) || { echo "FAIL vdc: stopped in $d (see $work/$d/vdc.log)"; fail=1; }
   done
   products=$(cd "$work/snap-full-1" && ls | grep -vxF -f ../snap-full-1.before | grep -vx vdc.log || true)
   bad=""
   for f in $products; do cmp -s "$work/snap-window-1/$f" "$work/snap-full-1/$f" || bad="$bad $f"; done
   if [ -n "$products" ] && [ -z "$bad" ]; then
      echo "PASS window snapshots through vdc: $(echo $products | wc -w) files identical"
   else
      echo "FAIL window snapshots through vdc: ${bad:- no files written}"; fail=1
   fi
else
   echo "vdc not found ($vdc): skipped the vdc comparison"
fi

# 4. run times
for m in off window; do
   echo "$m: $(cat "$work"/$m-[0-9]*/seconds | awk '{t += $1} END {printf "%.1f s per run", t/NR}')"
done

# 5. speed (PERF=<seconds> per run): a 300 x 20 x 20 nm film (FePt 10 nm high,
#    sensor layer above, laser FWHM 30 nm), written for PERF seconds with and
#    without the window; wall time per bit from the log, from bit 4 on (all
#    grains are integrated for the first 10 ps)
if [ "${PERF:-0}" != 0 ]; then
   sed -e 's/^material\[4\]:density = 1.0/material[4]:density = 0.01/' \
       -e 's/initial-spin-direction = 0.3,0,1/initial-spin-direction = 0,0,1/' \
       "$here/cases/hamr-granular/FePt_sensor.mat" > "$work/perf.mat"
   for m in off window; do
      mkdir "$work/perf-$m"
      cp "$work/perf.mat" "$work/perf-$m/FePt_sensor.mat"
      sed -e 's/system-size-x = 8 !nm/system-size-x = 300 !nm/' \
          -e 's/system-size-y = 8 !nm/system-size-y = 20 !nm/' \
          -e 's/laser-FWHM-\([xy]\) = 5/laser-FWHM-\1 = 30/' \
          -e 's/head-speed = 5000/head-speed = 1000/' \
          -e 's/head-field-x = 5/head-field-x = 14/' \
          -e 's/head-field-y = 6/head-field-y = 18/' \
          -e 's/bit-size = 4 !nm/bit-size = 6 !nm/' \
          -e 's/track-size = 7 !nm/track-size = 18 !nm/' \
          -e 's/number-of-bits = 1$/number-of-bits = 49/' \
          -e 's/output-rate = 500/output-rate = 1000/' \
          "$here/cases/hamr-granular/input" > "$work/perf-$m/input"
      [ $m = window ] && echo "hamr:active-window" >> "$work/perf-$m/input"
      (cd "$work/perf-$m" && timeout "$PERF" $mpi "$exe" > stdout 2>&1) || true
   done
   python3 - "$work" << 'EOF'
import sys, re, datetime
for m in ('off', 'window'):
   log = open(f'{sys.argv[1]}/perf-{m}/log').read()
   atoms = re.search(r'Number of atoms generated: (\d+)', log)
   bits = [(int(b), datetime.datetime.strptime(t, '%d-%m-%Y %H:%M:%S')) for t, b in re.findall(r'(\S+ \S+)\] Written bit (\d+) ', log.replace('[', ''))]
   bits = [b for b in bits if b[0] >= 4]
   window = [int(a) for a in re.findall(r'HAMR active window: \d+ grains, (\d+) atoms', log)]
   rate = f'{(bits[-1][1]-bits[0][1]).total_seconds()/(bits[-1][0]-bits[0][0]):.1f} s per bit (bits {bits[0][0]}-{bits[-1][0]})' if len(bits) > 1 else 'too few bits written, increase PERF'
   extra = f', integrated atoms {min(window)}-{max(window)}' if window else ''
   print(f'speed {m}: {rate}, {atoms.group(1) if atoms else "?"} atoms{extra}')
EOF
fi

exit $fail
