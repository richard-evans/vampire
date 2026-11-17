# VAMPIRE validation

Checks that a vampire binary (serial, CUDA DP/SP, multi-GPU, ...) produces the correct results by comparing against
analytic solutions and committed reference outputs from the serial CPU code.

```
python3 developers/test/validation/vtest.py run --exe ./vampire-cuda --label gpu-dp [--precision sp] [--jobs N]
python3 developers/test/validation/vtest.py compare --test gpu-dp --ref cpu
```

- `run` copies each `cases/<case>/` to `results/<label>/<case>/<run>/`, appends `--extra` lines (and
  `sim:integrator-random-seed` for seeded cases) to the input, and runs the binary there.
- `compare` exits non-zero on any failure; passing cases are summarised in one line (`--verbose` shows every quantity). Labels are looked up in `results/` first, then `reference/`.
  `--paired RTOL` instead compares every run with the same run (same seed) of the reference row by row, for result sets
  with the same dynamics whose output differs only by rounding.
- `developers/tests.sh -v` runs `vampire-serial` and compares against `reference/cpu`.

## GPU validation

Which check to run:

| change | check |
| ------ | ----- |
| meant to leave every result unchanged (refactoring, optimisation with the same arithmetic, new options that are off by default) | `identity_compare.sh` against the last commit that passed the full checks; if every run is identical, the results of those checks carry over |
| touches the multi-GPU version (`src/mgpu`, halo, decomposition) | `mgpu_check.sh` as well |
| changes results (new physics, different arithmetic or random numbers) | `run_gpu_suite.sh`, `mgpu_check.sh`, and `ab_compare.sh` for statistically equivalent changes |
| anything that could change what an ordinary user of the serial or MPI version sees (input, output files, checkpoints, vdc) | `compat_check.sh` against the reference version (CPU only, no GPU needed) |
| any change to the CPU code | upstream's regression suite, `developers/tests.sh -r` and `-m 2` (about 30 s and 50 s); built with GCC 15 the serial cases `config/binary-continuous`, `integrate/constrained-monte-carlo`, `program/hamr-cool`, `vio/all-columns`, `vio/grain-columns` and the 2-process cases `config/binary-continuous`, `config/output-mode-mpi-io`, `program/hamr-cool`, `vio/grain-columns` fail for the reference version too (low order bits of their golden files), so only other failures point to a change |
| the HAMR active window (`src/hamr/active_window.cpp`, window snapshots) | `hamr_window_check.sh` (below) |
| new feature or validation case | all of the above once, before committing |

The GPU scripts start several validation runs at a time on the GPU (`GPU_JOBS`, default 4 processes, so 2 runs on 2
ranks or 1 run on 4 ranks); results do not depend on this, and benchmarks always run alone. Use `GPU_JOBS=1` if the GPU
runs out of memory.

`developers/test/validation/run_gpu_suite.sh [label-prefix]` builds the CUDA version in double and single precision, runs the
suite with each build (and the `hamr-perf` benchmark), and packs results, build logs and system information into
`developers/test/validation/results/gpu-suite-<host>-<date>.tar.gz`. Each build is also run with statistics on the CPU
(two seeds), which must agree with the same runs of the default build to 1e-10 (the dynamics are the same; only the
summation order of the statistics differs, giving differences below 1e-12), and the CUDA binary with
`gpu:acceleration=off` must reproduce the serial references byte for byte (one seed per case, without the long
granular cases). It then runs the HAMR chain of `hamr_chain.sh` (granular medium of `hamr-granular-t0` with grain rounding
0.5): the CPU path relaxes it by Monte Carlo and writes it from the relaxed state (`sim:load-checkpoint=restart`) as a
reference; each build writes it from the same relaxed state, which must agree with the reference within the
deterministic tolerance, and the CPU path loads the checkpoint of each GPU write in a zero step read-back simulation
(hierarchical dipole field, binary `file-per-node` configuration output), whose spins must be bit for bit those of the
write. It runs `make clean` before each build and at the end.

`developers/test/validation/profile_gpu.sh [case] [dp|sp]` profiles one case (default `hamr-perf`) with Nsight Systems.

`developers/test/validation/ab_compare.sh [ref-commit] [dp|sp] [seeds]` compares the CUDA build of a reference commit
(default `HEAD`) with the working tree, for changes expected to give statistically equivalent rather than identical
results (for example a different random number sequence): it runs every finite temperature LLG case with both builds
using more seeds (default 32), compares them with each other and with the CPU references, and times the
`hamr-perf` and `thermal-perf` benchmarks with both builds.

`developers/test/validation/identity_compare.sh [ref-commit]` (default `develop`) is for changes that must not change any
result: it builds the CUDA version of a reference commit and of the working tree in both precisions, runs every
case with both (once as committed and once with extra statistics, configuration and checkpoint output) and compares
every output file byte for byte, ignoring only comment lines and timings. It also runs the cases with
non-magnetic atoms with and without `gpu:skip-non-magnetic-atoms` (see below) and checks that the output and grain
files and the final spins of every magnetic atom are identical, and prints the run times of both.

`gpu:skip-non-magnetic-atoms` (GPU versions, off by default) leaves the atoms of `non-magnetic=keep` materials
unintegrated: they keep their initial spin directions, which only appear in configuration and checkpoint files. Such
atoms enter no field or statistic of the magnetic atoms (statistics and dipole fields leave them out), so all other
results are bit for bit identical; the single GPU version still draws their random numbers so that the other atoms
receive the same ones. The option is refused with energy statistics (which include these atoms), with vectorial or
tensorial exchange, and when any exchange interaction involves an atom of a `non-magnetic=keep` material.

## Multi-GPU validation

`run --ranks N` starts each run with `mpirun -np N` (the `MPIRUN` environment variable replaces `mpirun`, e.g.
`MPIRUN="mpirun --oversubscribe --bind-to none"`; without `--bind-to none` concurrent runs bind their ranks to the
same cores). `vtest.py configs --test A --ref B` compares the final spin configurations of two
result sets atom by atom (atoms are paired by their coordinates, since MPI runs list atoms in rank order) and reports
whether they are bit for bit identical; runs need `--extra config:atoms=end --extra config:output-format=binary`.

**Bit for bit identical results on different numbers of GPUs require `gpu:decomposition-independent-order`**
(multi-GPU version only, off by default). It sums the exchange contributions of every atom in an order given by the
atom's position in the crystal (its unit cell and its number within the unit cell) instead of its local atom number,
which depends on the decomposition. Without it, runs on different numbers of GPUs agree only to floating point
round-off. It costs a sort of the neighbour list at start-up and has no measurable effect on the time per step.

`developers/test/validation/mgpu_check.sh [ranks ...]` (default `1 2 4`) builds the CPU MPI version and the multi-GPU version
(`make cuda-multigpu`) in both precisions and runs the zero temperature cases on each number of ranks, which may share
one GPU. One-rank runs are compared with the serial references and multi-rank runs with the CPU MPI version on the
same number of ranks (most cases start from random spins drawn per rank). `t0-spin-spiral` starts from a state that
does not depend on the decomposition: it is compared with the serial references on every number of ranks, and its
final spins must be bit for bit identical on every number of ranks (the script enables
`gpu:decomposition-independent-order` for the GPU runs). `t0-exchange-tensor` also starts from that state, with
tensorial exchange whose order of summation depends on the decomposition, so it is only compared with the references. The other LLG cases whose initial state does not depend on
the decomposition (finite temperature, HAMR and program cases) are run with one seed on every number of ranks and
their final spins must also be bit for bit identical, since the multi-GPU version draws the thermal noise of each atom
from its key and the time step; with all their seeds on one rank they are compared statistically with the serial
references. Finally `timeseries-T300` and `hamr-small` are rerun on every number of ranks saving a single checkpoint part way
(`sim:save-checkpoint=continuous` with a rate of 60 % of the time steps) and continued from it
(`sim:load-checkpoint=continue`); the final spins of the continuation must be bit for bit those of the uninterrupted
run. The granular HAMR writes (`hamr-granular`, `hamr-granular-t0`) save a checkpoint at the end, from which a one step
300 K time series is restarted (`sim:load-checkpoint=restart`, binary `file-per-node` configuration output, as in a
read-back simulation) on every number of ranks; its final spins must be bit for bit identical on every number of ranks,
and `hamr-granular-t0` (including its grain statistics) must agree with the first number of ranks within the
deterministic tolerance. Finally the granular writes are repeated with `gpu:skip-non-magnetic-atoms` on every number of
ranks: output and grain files and the final spins of every magnetic atom must be identical to those without it.
The HAMR chain of `hamr_chain.sh` (see `run_gpu_suite.sh` above) is run on every number of ranks, with the CPU MPI
version for the relaxation, the reference write and the loads of the multi-GPU checkpoints.
The checkpoints saved part way by `timeseries-T300` and `hamr-small` are also continued on another number of ranks
(each number of ranks loads those of the previous one in the list, cyclically): the atoms are found by their keys and
the final spins must still be bit for bit those of the uninterrupted run. Last, the chain is run on a medium with a
1 % dense fill between the grains (`material:density = 0.01` of the fill material, as in production runs) with
`create:decomposition-independent-random`, each stage on a different number of ranks: relaxation on n, multi-GPU
write on the next number of ranks (compared with the CPU write on n), CPU load of its checkpoint on the one after.

A checkpoint (one file per process) can be loaded on a different number of processes: the loader identifies atoms by
their key and stops with an error unless every atom is found exactly once with the same material. The structure must
therefore not depend on the decomposition. Random alloys, dilution (including the fill of `fill-space` materials) and
intermixing draw their random numbers per process unless `create:decomposition-independent-random` (off by default)
draws them from the seed and the atom's key instead. With `sim:load-checkpoint=continue` on a different number of
processes, the CPU versions restart their random number generator from the seed and the time step (statistically
equivalent but not identical to an uninterrupted run); the multi-GPU version draws the thermal noise from the atom key
and the time step, so with `gpu:decomposition-independent-order` its continuation is bit for bit that of the
uninterrupted run.

`run --seeds N` overrides the number of seeds of statistical cases. `compare --summary` prints the signed z-scores
of all statistical comparisons pooled together: when test and reference sample the same distribution they are
approximately standard normal (mean ~0, rms ~1, ~4.6 % beyond 2), which is a sharper test than the per-quantity
pass/fail. Rows of a per-row case are correlated, so judge the pooled values rather than a single case.

## Dipole solvers

`developers/test/validation/dipole_check.sh [vampire binary] [MPI ranks]` (CPU, a few seconds) runs a uniformly
magnetised Co film (40 x 12 x 3 nm, 1 nm cells) with the tensor solver (with and without
`dipole:on-the-fly-tensors`, which must give identical fields) and the hierarchical solver: and the hierarchical field must agree with the
tensor field to `TOL` tesla (default 1e-3; the film is longer than the range the hierarchical solver sums cell by cell,
so its far-field approximation is tested).

## HAMR active window

`developers/test/validation/hamr_window_check.sh [vampire binary]` checks `hamr:active-window` on small granular HAMR
films: with every grain in the window the output is identical to the run without it; over `SEEDS` seeds (default 16)
the magnetisation with and without the window agrees statistically (z-scores below `ZMAX`, default 4); window spin
snapshots rebuilt into full ones equal the full snapshots, and `vdc` (`VDC`, default `vdc/vdc` next to the binary)
writes identical files from both; and the run times. `FILM=small` runs every stage on a 30 x 8 nm film, small enough
for one CPU; `RANKS=N` runs the MPI and multi-GPU versions on N processes, `JOBS` runs several seeds at a time.

## Compatibility with a reference version

`developers/test/validation/compat_check.sh [ref-commit]` (default `upstream/v8`) builds the serial and MPI versions and
`vdc` of a reference commit and of the working tree, and checks on the CPU that the working tree behaves like the
reference for an ordinary user: every output file must be identical byte for byte (comment lines starting with `#`,
and the stdout and log files, are ignored). Its stages (`ONLY=stage,stage` selects some of them):

- `samples`: every sample of `samples/` with both serial versions (a sample that needs a set-up step, such as the
  unit cell generator, gets it; a run longer than 50000 steps is shortened for both versions)
- `config`: a small system with two magnetic and one non-magnetic material, written with every combination of
  `config:output-format` and `config:output-mode` (and one and two output nodes), serial and on `RANKS` MPI ranks
  (default 2); combinations that both versions refuse are reported. Every file, and every product of each version's own
  `vdc` (no option, `--xyz`, `--vtk`, `--povray`, `--text`, `--cells`, `--ssc`, `--grains`, `--povray-cells`,
  `--povray-sticks`, `--cell-size`) must be identical; for legacy files, which `vdc` cannot read, both must give the
  same error
- `old-files`: `vdc` of the working tree on the files written by the reference version must give the products of the
  reference `vdc`
- `checkpoints`: a 300 K run writes a checkpoint half way; it is continued (`sim:load-checkpoint=continue`) and
  restarted (`restart`) by every combination of writing and reading version, serial and parallel, and must give the
  results of the reference version reading its own checkpoint
- `window-snapshots` (working tree only): a granular HAMR write with `hamr:active-window` and spin configuration
  files, with and without `hamr:active-window-full-snapshots`: output and the `vdc` products (`--xyz`, `--vtk`,
  `--text`, `--povray`) must be identical
- `warnings`: the compiler warnings of each target; fails if the working tree warns in a file the reference does not

The binaries of a reference commit are cached in `results/compat-cache` (`CACHE`) and reused. Other variables: `JOBS`
(build jobs), `MPIRUN`, `WORK` (scratch directory, kept; by default a temporary directory, kept with `KEEP=1`). The
summary of PASS and FAIL lines and the times is written to `results/compat-<ref>-<date>/summary.txt`.

## Case format (`cases/<case>/case.json`)

All case inputs use `output:column-headers` and `output:precision=16`

| kind | runs | pass criterion |
| ------ | ------ | ---------------- |
| `analytic` / `precession` | 1 | max abs error vs analytic damped precession <= `tolerance[precision].atol` |
| `analytic` / `thermal` | `seeds` | seed-mean of time-averaged <mz> (after `skip_rows`) within `k`*stderr + `atol` of the Boltzmann average |
| `deterministic` (T=0) | 1 | every row: abs(test-ref) <= abs + atol*max(abs(ref)) + rtol*abs(ref) (`tolerance[precision]`, `abs` optional, default 0, for quantities close to zero in every compared row); an optional `rows` entry limits the comparison to the first rows, used in single precision for random initial states where round-off is amplified chaotically; reports "identical" when bit-equal |
| `statistical` (T>0) | `seeds` | abs(mean_test - mean_ref) <= `k` * sqrt(se_test^2 + se_ref^2) over seed time-averages (after `skip_rows`); with `"per_row": true` every output row is compared separately (e.g. one row per temperature) |
| `benchmark` | 1 | not compared, reports simulation run time; only run when named with `--cases` (e.g. `hamr-perf` for profiling) |

`"grain": true` adds the quantity `grain`: every column of `grain.txt` after time and temperature (grain and
material-grain statistics), compared element by element like the other quantities. References keep `grain.txt` for
these cases. `vtest.py configs` reads both the single file and the `file-per-node` configuration formats; with
`--magnetic-only` it ignores the atoms of the case's `non-magnetic=keep` materials. `vtest.py identical --test A --ref B`
compares the `output` and `grain.txt` files of every run of two result sets byte for byte, ignoring comment lines.

Precision (`dp`/`sp`) is taken from the test run's `meta.json` (set with `run --precision`).
Statistical checks fail spuriously ~0.1 % of the time per quantity: re-run with `--seed-offset 100` before treating a failure as real.

## Cases

- `analytic-*`: physics vs exact results (single-spin precession; Boltzmann average for non-interacting spins).
  `analytic-thermal-mc-anis` adds fourth and sixth order uniaxial anisotropy and temperature rescaling (Monte Carlo).
- `t0-*`: deterministic T=0 trajectories. `t0-anis-<term>` isolates one anisotropy term each; `biaxial`, `triaxial`
  and `cubic-rotated` are not implemented in the CUDA module, so they are expected to fail on GPU until fixed.
- `timeseries-T300`, `curie-llg`, `curie-mc`: thermal statistics (LLG and Monte Carlo are each compared only with
  the same integrator).
- `hamr-small` (statistical), `hamr-t0` (deterministic), `hamr-perf` (benchmark, ~168k atoms); `thermal-perf` (benchmark, the same system as a 300 K time series).
- `hamr-granular-t0` (deterministic, with grain statistics) and `hamr-granular` (statistical): 61,440 atom granular
  FePt film (rocksalt supercell, voronoi grains with 1 nm spacing, periodic in x and y, 8x8x20 nm, 27,160 FePt atoms)
  with a non-magnetic fill between the grains and a dense non-magnetic sensor layer above them (`non-magnetic=keep`),
  written by a HAMR head at 0 K and with a 300-800 K laser. Initial spins are tilted and the fill has full
  density, so that neither depends on the decomposition.
