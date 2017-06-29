#--------------------------------------------------------------
#                 Makefile for create module
#--------------------------------------------------------------

# List module object filenames
create_objects=\
create_system2.o \
cs_create_crystal_structure2.o \
cs_create_system_type2.o \
cs_create_neighbour_list2.o \
cs_set_atom_vars2.o \
alloy.o \
bubble.o \
bulk.o \
cone.o \
compare_radius.o \
cube.o \
cylinder.o \
data.o \
ellipsoid.o \
faceted.o \
initialize.o \
interface.o \
layers.o \
multilayers.o \
roughness.o \
sphere.o \
teardrop.o \
truncated_octahedron.o \
voronoi.o \
voronoi_grain_rounding.o \
voronoi_substructure.o \
voronoi_vertex_points.o

# Append module objects to global tree
OBJECTS+=$(addprefix obj/create/,$(create_objects))
