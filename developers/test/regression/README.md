# Regression tests

Golden-output tests that detect *any* change in results, for the serial and MPI builds separately.
They do not check that the results are physically right; that is the job of the integration tests.

```bash
./run_regression.sh                      # vampire-serial, all cases
./run_regression.sh --mpi 2              # vampire-parallel on 2 processes
./run_regression.sh --case 'integrate/*' # only cases matching a glob
./run_regression.sh --exact              # byte-for-byte, for changes that must not alter any bits
```
or from the repository root, `developers/tests.sh -r` (serial) and `developers/tests.sh -m [N]` (MPI).

## Cases

Each case is a directory `cases/<module>/<case>/` containing:

| file | content |
|---|---|
| `input`, `*.mat`, ... | inputs, copied to a temporary directory for the run |
| `expected-serial.dat` | golden `output` file from `vampire-serial` (comment lines removed) |
| `expected-mpi-np<N>.dat` | golden `output` file from `vampire-parallel` on N processes |
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
