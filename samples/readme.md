# Vampire samples

Example simulations, roughly in the order a new user might work through them.
Each folder has an `input` file, the material files it uses and a `readme.md`.
Run a sample from inside its folder, for example

```bash
cd samples/geometry
../../vampire-serial
```

Output (`output`, `log` and any configuration files) is written to the same folder.
Configuration files can be converted for visualisation with the vampire data
converter, e.g. `../../vdc/vdc --povray` (see the manual for the options).

To check that vampire has compiled correctly, run the `input` and `Co.mat` in the
top-level folder first (a short benchmark, a few seconds).

| Sample | Shows | Run time |
|--------|-------|----------|
| [spin-textures/domain-wall](spin-textures/domain-wall) | initial spin configuration with a 180 degree domain wall | < 1 s |
| [spin-textures/skyrmion](spin-textures/skyrmion) | initial spin configuration with a skyrmion | < 1 s |
| [spin-textures/vector-field](spin-textures/vector-field) | initial spin configuration interpolated from a user vector field | < 1 s |
| [geometry](geometry) | cutting materials to 2D shapes with geometry files | < 1 s |
| [grain-alternating-magnetisation](grain-alternating-magnetisation) | Voronoi granular film with alternating grain magnetisation | ~1 s |
| [unitcells](unitcells) | creating a custom unit cell file (Python or C++) for a ferrimagnet | ~6 s |
| [laser-electrical-pulse](laser-electrical-pulse) | laser pulse (two-temperature model) followed by a spin-transfer torque pulse | < 1 s |
| [sld](sld) | spin-lattice dynamics in bcc Fe | ~40 s |
