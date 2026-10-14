# Geometry files

Shows `material:geometry-file`, which cuts a material to a 2D shape, similar to
lithography. A 20 x 20 x 2 nm film has two layers:

- material 1 (lower half, `minimum-height = 0.0`, `maximum-height = 0.5`) is cut to the triangle in `triangle.geo`
- material 2 (upper half, `minimum-height = 0.5`, `maximum-height = 1.0`) is cut to the rectangle in `rectangle.geo`

A geometry file gives the number of points on the first line, then one `x y`
point per line as a fraction (0-1) of the system size. The last point joins up
with the first.

```bash
cd samples/geometry
../../vampire-serial
../../vdc/vdc --xyz
```

The log reports how much of the system each material makes up (about 21% triangle,
79% rectangle). `vdc --xyz` writes the structure to `crystal.xyz`, which can be
opened in a structure viewer such as Jmol or OVITO. The two layers are given
different elements (Co, Fe) so they show in different colours there.
