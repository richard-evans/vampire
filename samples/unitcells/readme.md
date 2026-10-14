# Custom unit cells

Vampire can read any crystal structure from a unit cell file (`.ucf`), which
lists the atoms in the unit cell and the exchange interactions between them.
This example creates one for a rare earth - transition metal (RE-TM)
ferrimagnet: an fcc lattice with one RE atom and three TM atoms per unit cell,
ferromagnetic RE-RE and TM-TM exchange and antiferromagnetic RE-TM exchange.

The unit cell creator finds all nearest neighbour interactions and writes
`UC.ucf` and a matching material file `MAT.mat`. There are two versions that
give identical files:

- `unit-cell-creator.py` (Python 3, no extra packages): `python3 unit-cell-creator.py`
- `unit-cell-creator.cpp`: `g++ -o unit-cell-creator unit-cell-creator.cpp && ./unit-cell-creator`

Then run vampire:

```bash
cd samples/unitcells
python3 unit-cell-creator.py
../../vampire-serial
```

`input` runs 10 ps of spin dynamics at 300 K for a 5 nm cube and outputs the
magnetisation of each material. The two sublattices stay antiparallel.

To make your own structure, edit the parameters at the top of the script: the
unit cell size, the materials, the exchange constants between materials, the
atoms in the unit cell (fractional coordinates and material number) and the
interaction range. The material file sets `material[n]:unit-cell-category = n`
so that atoms with material number n-1 in the unit cell file belong to material
n. The format of the unit cell file is described in the manual
(`manual/unit-cell-file`).
