# Spin textures

Three examples of `material:initial-spin-direction`, which sets up a
non-uniform starting spin configuration in a 20 x 20 x 1 nm cobalt film.
The texture is set in each folder's `Co.mat`; the `input` files are the same.

| Folder | Setting in Co.mat |
|--------|-------------------|
| `domain-wall` | `initial-spin-direction = domain-wall-x, 0.5, 0.2` (centre, width as fractions of the system size) |
| `skyrmion` | `initial-spin-direction = skyrmion, 0.5, 0.5, 0.3` (centre x, y and radius) |
| `vector-field` | `initial-spin-direction = vector-field, field.dat` (points and directions in `field.dat`) |

Each runs a single time step and writes the spin configuration (`config:atoms`).
To see the texture:

```bash
cd samples/spin-textures/skyrmion
../../../vampire-serial
../../../vdc/vdc --povray
```

then render `spins.pov` with POV-Ray. Change the time steps and temperature in
`input` to watch the texture evolve.
