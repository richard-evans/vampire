Vampire
==============

By Richard Evans

Vampire is a high performance general purpose code for the atomistic simulation of magnetic materials. Using a variety of common simulation methods it can calculate the equilibrium and dynamic magnetic properties of a wide variety of magnetic materials and phenomena, including ferro, ferri and antiferromagnets, core-shell nanoparticles, ultrafast spin dynamics, magnetic recording media, heat assisted magnetic recording, exchange bias, magnetic multilayer films and complete devices.

Getting started
---------------
Vampire needs a C++17 compiler; MPI is only needed for the parallel version.

**Build**
```bash
make serial          # vampire-serial, using g++
make serial-llvm     # vampire-serial, using clang++ (recommended on macOS, where g++ is clang)
make parallel        # vampire-parallel, using mpicxx (needs MPI)
make vdc             # vdc, the vampire data converter for visualisation
```
`make` on its own builds all of serial, parallel and vdc. Add `-j 8` to compile in parallel.

**Check the build**

The `input` and `Co.mat` files in this folder run a short benchmark (a few seconds):
```bash
./vampire-serial
```
It prints the magnetisation as it runs and writes the results to `output` and a `log` file.
The parallel version runs the same way under MPI, e.g. `mpirun -np 4 ./vampire-parallel`.

**Run the samples**

The `samples` folder has example simulations with a readme for each: spin textures, geometry,
granular films, custom unit cells, laser and spin-transfer torque pulses and spin-lattice dynamics.
Start with `samples/readme.md`.

**Visualise**

With `config:atoms` set in the input file, vampire writes the spin configurations, which `vdc`
converts for POV-Ray, VTK (ParaView) and structure viewers, e.g. `vdc/vdc --povray`.

**Documentation**

The full documentation, including every input and material file keyword, is on the vampire website
https://vampire.york.ac.uk and in the LaTeX manual in `manual/`.

**Developers**

Tests, developer tools and notes for working on the code are in `developers/`. Run the test
suite with `developers/tests.sh`.

Capabilities
---------------
Vampire is designed to be highly flexible to deal with a wide variety of problems using a diverse set of simulation tools and methods. The capabilities of the code can be summarised broadly in terms of the simulation methods, standard problems, structural properties and features of the code, all of which can be combined to tackle almost any problem.

**Simulation methods**
-Stochastic Landau-Lifshitz-Gilbert equation (spin dynamics)
-Monte Carlo metropolis
-Constrained Monte Carlo metropolis

**Standard calculations**
-Ultrafast spin dynamics
-Hysteresis loops
-Curie temperature
-Temperature dependent anisotropy
-Temperature dependent energy barriers
-Field cooling
-Heat assisted and conventional magnetic recording
-Laser induced spin dynamics

**Structural properties**
-Bulk-like systems
-Thin films
-Nanoparticles - spheres, cubes, truncated octahedra, cylinders
-Voronoi granular structures
-Nanoparticle arrays
-Core-shell nanoparticles
-Multilayer thin films
-Interface roughness and intermixing
-Dilute magnetic systems
-Lithographically defined geometries
-SC, FCC, HCP, and BCC crystal structures
-User-defined atomic structures - for example from Molecular Dynamics simulations

**Magnetic properties**
-Ferromagnets
-Antiferromagnets
-Ferrimagnets
-Spin glass
-Single-ion, 2-ion and cubic anisotropies
-Scalar, vector and tensor forms of exchange including the DM interaction
-User-defined Hamiltonian from ab-initio Density Functional Theory (DFT) calculations
-Demagnetisation fields (macrocell approximation)

**Code features**
-Modular object-oriented C++
-Simple to use textfile input
-High performance code
-Parallelisation using the MPI library
-Variety of geometric decomposition algorithms
-Usable on a laptop to a supercomputer with thousands of cores
-Output to PoVRAY for visualisation and publication quality graphics
-Output to rasmol/jmol for structural inspection
-Minimal dependence on external libraries for portability
-Freely available open source code

License
---------------
Vampire as a whole is distributed under the GNU GPL version 2 (see `licence-gpl`).
Most source files are individually available under the Free BSD licence (see `licence-bsd`);
the header of each file states which licence applies to it. Files derived from JAMS
(https://github.com/stonerlab/jams) also carry the original MIT licence notice in their header.

