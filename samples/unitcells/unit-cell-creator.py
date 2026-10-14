#!/usr/bin/env python3
#-------------------------------------------------------------------------------
#
#   Unit cell creator for vampire
#
#   Takes the atoms in a unit cell, finds all nearest neighbour exchange
#   interactions (including those with atoms in the neighbouring unit cells)
#   and writes a unit cell file (UC.ucf) and a matching material file (MAT.mat)
#   for vampire.
#
#   The example is a rare earth - transition metal (RE-TM) ferrimagnet on an
#   fcc lattice: one RE atom and three TM atoms per unit cell, with
#   antiferromagnetic RE-TM exchange. Edit the parameters below for your own
#   system and run with
#
#      python3 unit-cell-creator.py
#
#   Python version of unit-cell-creator.cpp (R F L Evans 2015).
#
#-------------------------------------------------------------------------------

import itertools

#-------------------------------------------------------------------------------
# Parameters
#-------------------------------------------------------------------------------

# unit cell size (Angstroms) in x, y and z
unit_cell_size = (3.54, 3.54, 3.54)

# materials: name, moment (Bohr magnetons), damping, uniaxial anisotropy
# (J/atom), element (used for visualisation) and initial spin direction
materials = [
   dict(name="RE", mu_s=7.63, alpha=0.1, Ku=-8.07246e-24, element="Ag", spin=(0.0, 0.0, -1.0)),
   dict(name="TM", mu_s=1.92, alpha=0.1, Ku=-8.07246e-24, element="Fe", spin=(0.0, 0.0,  1.0)),
]

# exchange constants (J/link) between materials, exchange[i][j] for i-j pairs
exchange = [
   [ 1.26e-21, -1.09e-21],  # RE-RE, RE-TM
   [-1.09e-21,  2.835e-21], # TM-RE, TM-TM
]

# atoms in the unit cell: fractional coordinates (cx, cy, cz), material,
# lattice category (lc) and height category (hc). Materials are numbered
# from 0 and index the materials list above.
unit_cell = [
   dict(cx=0.0, cy=0.0, cz=0.0, material=0, lc=1, hc=0),
   dict(cx=0.5, cy=0.5, cz=0.0, material=1, lc=0, hc=0),
   dict(cx=0.5, cy=0.0, cz=0.5, material=1, lc=0, hc=1),
   dict(cx=0.0, cy=0.5, cz=0.5, material=1, lc=0, hc=1),
]

# interaction range (in units of the unit cell size): atoms closer than this
# are neighbours. sqrt(0.5) is the nearest neighbour distance in fcc.
interaction_range = (0.5*0.5 + 0.5*0.5)**0.5

#-------------------------------------------------------------------------------
# Find all neighbours of the atoms in the unit cell, looking at atoms in the
# same unit cell and in the 26 unit cells around it (dx, dy, dz = -1, 0, 1)
#-------------------------------------------------------------------------------
interactions = []
for i, atom_i in enumerate(unit_cell):
   for dx, dy, dz in itertools.product((-1, 0, 1), repeat=3):
      for j, atom_j in enumerate(unit_cell):
         # skip the atom itself
         if i == j and (dx, dy, dz) == (0, 0, 0):
            continue
         rx = atom_j["cx"] + dx - atom_i["cx"]
         ry = atom_j["cy"] + dy - atom_i["cy"]
         rz = atom_j["cz"] + dz - atom_i["cz"]
         if rx*rx + ry*ry + rz*rz <= interaction_range**2:
            Jij = exchange[atom_i["material"]][atom_j["material"]]
            interactions.append((i, j, dx, dy, dz, Jij))

# vampire expects numbers in a compact general format, e.g. 0.5, 1, -1.09e-21
def fmt(x):
   return f"{x:g}"

#-------------------------------------------------------------------------------
# Write unit cell file
#-------------------------------------------------------------------------------
with open("UC.ucf", "w") as ucf:
   ucf.write("# Unit cell size:\n")
   ucf.write("\t".join(fmt(s) for s in unit_cell_size) + "\n")
   ucf.write("# Unit cell vectors: \n")
   ucf.write("1.0 0.0 0.0 \n")
   ucf.write("0.0 1.0 0.0 \n")
   ucf.write("0.0 0.0 1.0 \n")
   ucf.write("# Atoms num, id cx cy cz mat lc hc \n")
   ucf.write(f"{len(unit_cell)}\n")
   for n, atom in enumerate(unit_cell):
      ucf.write(f"{n}\t{fmt(atom['cx'])}\t{fmt(atom['cy'])}\t{fmt(atom['cz'])}\t"
                f"{atom['material']}\t{atom['lc']}\t{atom['hc']}\n")
   # number of interactions and exchange type
   ucf.write("#Interactions n exctype, id i j dx dy   dz        Jij\n")
   ucf.write(f"{len(interactions)}\tisotropic\n")
   for n, (i, j, dx, dy, dz, Jij) in enumerate(interactions):
      ucf.write(f"{n}\t{i}\t{j}\t{dx}\t{dy}\t{dz}\t{fmt(Jij)}\n")

#-------------------------------------------------------------------------------
# Write material file
#-------------------------------------------------------------------------------
line = "#------------------------------------------------\n"
with open("MAT.mat", "w") as mat:
   mat.write("#================================================\n")
   mat.write("# Generated material file for input into vampire\n")
   mat.write("#================================================\n")
   mat.write("#\n")
   mat.write("# File timestamp: \n")
   mat.write("#\n")
   mat.write(line)
   mat.write("# Number of Materials\n")
   mat.write(line)
   mat.write(f"material:num-materials={len(materials)}\n")
   mat.write(line)
   for m, material in enumerate(materials, start=1):
      sx, sy, sz = material["spin"]
      mat.write(f"# Material {m} ({material['name']})\n")
      mat.write(line)
      mat.write(f"material[{m}]:material-name={material['name']}\n")
      # atoms in the unit cell file with material m-1 belong to this material
      mat.write(f"material[{m}]:unit-cell-category={m}\n")
      mat.write(f"material[{m}]:damping-constant={fmt(material['alpha'])}\n")
      mat.write(f"material[{m}]:atomic-spin-moment={fmt(material['mu_s'])} !muB\n")
      mat.write(f"material[{m}]:uniaxial-anisotropy-constant={fmt(material['Ku'])}\n")
      mat.write(f"material[{m}]:material-element={material['element']}\n")
      mat.write(f"material[{m}]:initial-spin-direction={fmt(sx)},{fmt(sy)},{fmt(sz)}\n")
      mat.write(line)

print(f"Wrote UC.ucf ({len(unit_cell)} atoms, {len(interactions)} interactions) and MAT.mat ({len(materials)} materials)")
