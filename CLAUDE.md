# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Vampire is a C++ atomistic simulation code for magnetic materials. It solves the stochastic Landau-Lifshitz-Gilbert (LLG) equation and Monte Carlo methods to calculate equilibrium and dynamic magnetic properties. It supports MPI parallelism and CUDA acceleration. The OpenCL code is deprecated.

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
# Run all tests
./tests.sh -a

# Integration tests only
./tests.sh -i
# or: cd test/integration && make && ./integration_tests

# Unit tests only
./tests.sh -u
# or: cd test/unit && make && ./unit_tests
```

Integration tests are compiled separately with C++17 (`test/integration/makefile`). Unit tests link against a subset of Vampire object files (`test/unit/makefile`).

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
| `cuda` / `opencl` | GPU acceleration backends (compiled with `CUDA`/`OPENCL` defines) |
| `qvoronoi` | Embedded qhull library for Voronoi tessellation |

### Global Data

Atom data lives in flat SoA (Structure-of-Arrays) vectors in the `atoms` namespace (`hdr/atoms.hpp`):
- Coordinates: `atoms::x/y/z_coord_array`
- Spin vectors: `atoms::x/y/z_spin_array`
- Fields: `atoms::x/y/z_total_spin_field_array`, `x/y/z_total_external_field_array`, `x/y/z_thermal_field_array`
- Neighbour list: CSR format via `neighbour_list_array`, `neighbour_list_start_index`, `neighbour_list_end_index`
- Material/grain/cell indices: `type_array`, `grain_array`, `cell_array`

Simulation state is in the `sim` namespace (`hdr/sim.hpp`): integrator enum, program enum, `sim::time`, `sim::total_time`, `sim::partial_time`, temperature, field, etc.

Material parameters are in the `mp` namespace (`src/main/initialise_variables.cpp`), loaded from `.mat` files.

### Input Files

- `input` — main parameter file, keyword:value pairs with optional `!units` suffix (e.g. `dimensions:system-size-x = 7.7 !nm`)
- `*.mat` — material parameter files referenced from `input` via `material:file = Co.mat`

### Compiler Defines

| Define | Effect |
|--------|--------|
| `MPICF` | Enable MPI parallel code paths |
| `CUDA` | Enable CUDA GPU backend |
| `OPENCL` | Enable OpenCL GPU backend |
| `FFT` | Enable FFTW for spin waves and FFT dipole |
| `COMP` | Compiler identification string printed at startup |

### Adding a New Module

Use `util/initialise_new_vampire_module.cpp` as a template. Each module needs: `src/<name>/`, `hdr/<name>.hpp`, `src/<name>/internal.hpp`, `src/<name>/makefile`, and an entry in the root makefile's `include` list.
