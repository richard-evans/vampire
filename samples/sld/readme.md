# Spin-lattice dynamics

Coupled spin and lattice dynamics (`sim:integrator = spin-lattice`) in a
2.87 nm cube of bcc Fe (2000 atoms) at 300 K, using a harmonic interatomic
potential and pseudodipolar spin-lattice coupling. The spin-lattice parameters
(mass, lattice damping, exchange J0, potential V0, coupling C0) are in `fe.mat`.

```bash
cd samples/sld
../../vampire-serial
```

Takes about 40 s. `output` has columns for time, magnetisation, spin temperature
and lattice temperature. Both start near zero: the spin temperature rises
quickly to about 300 K, while the lattice heats up more slowly through the
spin-lattice coupling (about 140 K after 1 ps).

To continue a run, uncomment `sim:save-checkpoint=end` for the first run and
`sim:load-checkpoint=continue` for the next one.
