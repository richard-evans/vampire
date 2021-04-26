#--------------------------------------------------------------
#                 Makefile for create module
#--------------------------------------------------------------

# List module object filenames
create_objects=\
create.o \
cs_create_crystal_structure2.o \
cs_create_system_type2.o \
cs_set_atom_vars2.o \
alloy.o \
bubble.o \
bulk.o \
centre_particle.o \
cone.o \
compare_radius.o \
composition.o \
cube.o \
cylinder.o \
data.o \
ellipse.o \
ellipsoid.o \
faceted.o \
get_material_height.o \
hexagonal_array.o \
initialize.o \
interface.o \
layers.o \
multilayers.o \
mpi.o \
particle.o \
roughness.o \
sort_atoms_by_grain.o \
sphere.o \
square_array.o \
teardrop.o \
truncated_octahedron.o \
voronoi.o \
voronoi_grain_rounding.o \
voronoi_substructure.o \
voronoi_vertex_points.o

# Append module objects to global tree
OBJECTS+=$(addprefix obj/create/,$(create_objects))
