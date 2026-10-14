# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Vampire is a C++ atomistic simulation code for magnetic materials. It solves the stochastic Landau-Lifshitz-Gilbert (LLG) equation and Monte Carlo methods to calculate equilibrium and dynamic magnetic properties. It supports MPI parallelism and CUDA acceleration.

## Repository Layout

The top level is kept for users; developer-only material lives in `developers/`.

| Path | Contents |
|------|----------|
| `src/`, `hdr/`, `obj/` | Vampire source, public headers, object files |
| `vdc/` | Vampire data converter (separate program, `make vdc`) |
| `samples/` | Example simulations for users |
| `manual/` | LaTeX user manual |
| `developers/test/` | Unit, integration and regression tests |
| `developers/tests.sh` | Test runner |
| `developers/tools/` | Developer tools: module generator, style check, CUDA remote-run scripts, UCF utilities |

## Build Commands

```bash
# Build serial and parallel executables plus vdc data converter (default)
make

# Build only serial executable (vampire-serial)
make serial-llvm -j 8

# Build only MPI parallel executable (vampire-parallel)
make parallel-llvm -j 8

# Build the vampire data converter utility
make vdc

# Debug builds
make serial-debug       # GCC debug with profiling/coverage flags
make parallel-debug     # MPI debug build

# Clean object files only
make clean

# Remove object files and executables
make purge

# Install to /opt/vampire (requires sudo)
make install
```

Run the compiled binary with an input file in the working directory:
```bash
./vampire-serial
./vampire-serial --input-file my_input  # specify alternate input file name
./vampire-parallel                       # MPI build, run with mpirun
```

## Testing

```bash
# Run standard integration and unit tests (default, fast)
developers/tests.sh

# Run all tests, including extended (long running) integration tests
developers/tests.sh -a

# Standard integration tests only
developers/tests.sh -i
# or: cd developers/test/integration && make && ./integration_tests

# Extended integration tests only (statistical tests against exact results, several seconds each)
developers/tests.sh -e
# or: ./integration_tests --extended   (--all runs standard and extended)

# Unit tests only
developers/tests.sh -u
# or: cd developers/test/unit && make && ./unit_tests
```

Integration tests are compiled separately with C++17 (`developers/test/integration/makefile`). Unit tests link against a subset of Vampire object files (`developers/test/unit/makefile`).

## Architecture

### Execution Flow

```
main() → mp::initialise()  // parse input and material files
       → cs::create()      // build atomic structure (geometry, lattice, neighbours)
       → sld::initialize() // spin-lattice dynamics initialisation
       → sim::run()        // dispatch to program type, run integrator loop
```

### Module Layout (`src/` and `hdr/`)

Each module lives in `src/<name>/` with a corresponding `hdr/<name>.hpp` public header. Module-private declarations go in `src/<name>/internal.hpp`. The root `makefile` includes each module's `src/<name>/makefile`.

Key modules:

| Module | Purpose |
|--------|---------|
| `main` | Entry point, command-line parsing (`--input-file`, `--output-file`, `--version`), calls `mp::initialise` |
| `material` | Material container (`mp::materials_t`), `.mat` file reader, dispatch of material keywords to each module's `match_material_parameter`, core material keywords (damping, moment, name, etc.) |
| `create` | Atomic structure generation: lattice types, nanoparticle shapes, Voronoi granular structures, multilayers |
| `simulate` | Integrators (LLG-Heun, LLG-midpoint, Monte Carlo, CMC, LLB, LSF); dispatches to `program::` |
| `program` | High-level simulation programs (benchmark, hysteresis, Curie temperature, HAMR, field sweep, etc.) |
| `exchange` | Exchange interaction types: scalar, vector (DM), tensor, biquadratic, four-spin |
| `anisotropy` | Uniaxial, cubic, surface Néel anisotropy energies and fields |
| `dipole` | Macrocell demagnetisation field calculation |
| `statistics` | Online statistics (magnetisation, susceptibility, etc.) — called each `partial_time` |
| `vio` | All I/O: input file parsing (`vin::match`), output (`vout::data`), unit conversion, logging |
| `mpi` | MPI domain decomposition wrappers |
| `cells` | Macrocell partitioning used by dipole and statistics |
| `unitcell` | Unit cell definitions (SC, FCC, BCC, HCP) and UCF file reader |
| `neighbours` | Neighbour list construction (cutoff radius) |
| `ltmp` | Localised temperature model |
| `hamr` | Heat-assisted magnetic recording programs |
| `spintorque` | Spin transfer torque |
| `spinwaves` | Spin wave dispersion (requires FFTW: uncomment in makefile) |
| `micromagnetic` | Micromagnetic hybrid mode |
| `gpu` | Backend-independent wrapper that dispatches to a GPU acceleration backend (currently `cuda`) |
| `cuda` | CUDA GPU acceleration backend (compiled with the `CUDA` define) |
| `atoms` | Global atom data (`atoms::` SoA arrays declared in `hdr/atoms.hpp`) |
| `random` | Random number generation (`mtrandom::`, MT19937 in `mtrand.cpp`) |
| `utility` | Support library: errors (`err::vexit`, `err::fatal`, `err::warning`), unit conversion (`units::`), vector maths (`vmath::`) |

### Global Data

Atom data lives in flat SoA (Structure-of-Arrays) vectors in the `atoms` namespace (`hdr/atoms.hpp`):
- Coordinates: `atoms::x/y/z_coord_array`
- Spin vectors: `atoms::x/y/z_spin_array`
- Fields: `atoms::x/y/z_total_spin_field_array`, `x/y/z_total_external_field_array`, `x/y/z_thermal_field_array`
- Neighbour list: CSR format via `neighbour_list_array`, `neighbour_list_start_index`, `neighbour_list_end_index`
- Material/grain/cell indices: `type_array`, `grain_array`, `cell_array`

Simulation state is in the `sim` namespace (`hdr/sim.hpp`): integrator enum, program enum, `sim::time`, `sim::total_time`, `sim::partial_time`, temperature, field, etc.

Material parameters are in the `mp` namespace (`hdr/material.hpp`), loaded from `.mat` files by the `material` module into `mp::read_material` and then copied to `mp::material`. Global initialisation and derived parameters (`mp::initialise`, `mp::set_derived_parameters`) remain in `src/main/initialise_variables.cpp`.

### Input Files

- `input` — main parameter file, keyword:value pairs with optional `!units` suffix (e.g. `dimensions:system-size-x = 7.7 !nm`)
- `*.mat` — material parameter files referenced from `input` via `material:file = Co.mat`

### Compiler Defines

| Define | Effect |
|--------|--------|
| `MPICF` | Enable MPI parallel code paths |
| `CUDA` | Enable CUDA GPU backend |
| `FFT` | Enable FFTW for spin waves and FFT dipole |
| `COMP` | Compiler identification string printed at startup |

### Adding a New Module

Use `developers/tools/initialise_new_vampire_module.cpp` as a template. Each module needs: `src/<name>/`, `hdr/<name>.hpp`, `src/<name>/internal.hpp`, `src/<name>/makefile`, and an entry in the root makefile's `include` list.
