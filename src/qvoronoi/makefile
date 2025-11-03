#--------------------------------------------------------------
#          Makefile for qvoronoi module
#--------------------------------------------------------------

# List module object filenames
QV_OBJECTS=\
geom.o\
geom2.o\
global.o\
io.o\
libqhull.o\
mem.o\
merge.o\
poly.o\
poly2.o\
qhrandom.o\
qset.o\
qvoronoi.o\
rboxlib.o\
stat.o\
user.o\
usermem.o\
userprintf.o\
userprintf_rbox.o

# Append module objects to global tree
OBJECTS+=$(addprefix obj/qvoronoi/,$(QV_OBJECTS))

# Add specific options for qvoronoi to disable annoying warnings
obj/qvoronoi/%.o : OPTIONS = -w
