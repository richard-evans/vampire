# Regression tests

Golden-output tests that detect *any* change in results, for the serial and MPI builds separately.
They do not check that the results are physically right; that is the job of the integration tests.

```bash
./run_regression.sh                      # vampire-serial-test, all cases
./run_regression.sh --mpi 2              # vampire-parallel-test on 2 processes
./run_regression.sh --case 'integrate/*' # only cases matching a glob
./run_regression.sh --exact              # byte-for-byte, for changes that must not alter any bits
./run_regression.sh --exe ./vampire-serial  # another executable, eg an optimised build
```
or from the repository root, `developers/tests.sh -r` (serial) and `developers/tests.sh -m [N]` (MPI), which build the
test executables first.

## Test builds

The goldens are made with the reproducible test builds, `make serial-test` and `make parallel-test`
(`vampire-serial-test`, `vampire-parallel-test`). They turn off fused multiply-add contraction (`-ffp-contract=off`) and
use no architecture-specific flags, so they give the same bits with clang on arm64 and gcc on x86-64. The optimised
builds (`serial`, `serial-llvm`, ...) fuse `a*b+c` wherever the target has FMA instructions (clang on Apple silicon
always does). This rounds once instead of twice, which is no less accurate, but it changes the last bits, and thermal and
chaotic cases grow them into different trajectories. An optimised build may therefore fail some cases here; check its
physics with the validation suite (`developers/test/validation`), which uses tolerances and statistics.

## Cases

Each case is a directory `cases/<module>/<case>/` containing:

| file | content |
|---|---|
| `input`, `*.mat`, ... | inputs, copied to a temporary directory for the run |
| `expected-serial.dat` | golden `output` file from `vampire-serial-test` (comment lines removed) |
| `expected-mpi-np<N>.dat` | golden `output` file from `vampire-parallel-test` on N processes |
| `tolerance` | optional, `rel <r> abs <a>`; the default is `rel 1e-12 abs 1e-14` |
| `files` | optional, other output files to compare, one per line (glob patterns allowed, as serial and MPI runs can write different files), with goldens `expected-serial.<file>` and `expected-mpi-np<N>.<file>`; binary files are compared byte for byte |
| `expected-error` | for input that must be rejected: vampire must fail and its terminal output must contain each line of this file; runs in every mode |

Each run is stopped after `TIME_LIMIT` seconds (default 300), so a hang counts as a failure. A missing `output` file (some programs never write one) is compared as an empty file.
A case runs only in the modes for which it has a golden file. Numeric fields pass if
`|obtained - expected| <= abs + rel*|expected|`.

## Adding a case

1. Make `cases/<module>/<name>/` with a minimal `input` and `.mat` file: under about 1k atoms and 2 s,
   with `sim:integrator-random-seed` and `output:precision = 16` set.
2. Generate the golden files with `./run_regression.sh --bless --case <module>/<name>`, and with
   `--mpi 2` added if the case runs in parallel.
3. Check the golden output is physically sensible, then commit the case and its golden files together.

Only re-bless an existing case when a change in results is intended and understood, and say why in the
commit message.
