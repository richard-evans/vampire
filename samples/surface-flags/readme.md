# Surface flags in a unit cell file

Néel pair anisotropy acts on surface atoms. By default vampire finds them
geometrically: an atom with fewer neighbours than a threshold is a surface atom.
For a complete system read from a unit cell file this is not always what is
wanted, so the unit cell file can give an optional last column on each atom
line, 1 for a surface atom and 0 for a bulk atom (given for every atom or none).

`island.ucf` holds a 6 x 6 x 3 simple cubic island (108 atoms) as a single unit
cell. The top and side faces are flagged as surface atoms. The bottom face rests
on a substrate, so it is flagged as bulk even though its atoms are missing a
neighbour. Néel anisotropy then acts on the 76 flagged atoms; the geometric
rule with a threshold of six neighbours would give 92.

```bash
cd samples/surface-flags
../../vampire-serial
```

The log reports `Identifying surface atoms using unit cell surface flags.` and
`76 surface atoms found.` `input` runs 20 ps of spin dynamics at 10 K and outputs
the magnetisation and anisotropy energy. The unit cell file format, including the
surface flag column, is described in the manual (`manual/unit-cell-file`).
