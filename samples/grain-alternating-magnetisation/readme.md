# Granular film with alternating grain magnetisation

A 30 x 30 x 2 nm cobalt film is divided into Voronoi grains (`create:voronoi-film`,
5 nm grains with 1 nm spacing). With
`create:grain-magnetisation-direction = alternating`, odd-numbered grains start
with their magnetisation reversed, giving a chessboard-like pattern of up and
down grains.

```bash
cd samples/grain-alternating-magnetisation
../../vampire-serial
../../vdc/vdc --povray
```

Runs a single time step and writes the spin configuration for viewing.
