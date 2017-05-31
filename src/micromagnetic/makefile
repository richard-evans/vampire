#--------------------------------------------------------------
#          Makefile for micromagnetic module
#--------------------------------------------------------------

# List module object filenames
micromagnetic_objects =\
data.o \
initialize.o \
interface.o \
calculate_A.o \
calculate_alpha.o \
calculate_chi_para.o \
calculate_chi_perp.o \
calculate_gamma.o \
calculate_ku.o \
calculate_ms.o \
calculate_tc.o \
atomistic_llg_heun.o \
micromagnetic_llb.o \
micromagnetic_llg.o \
fieldsllb.o \
fieldsllg.o

# Append module objects to global tree
OBJECTS+=$(addprefix obj/micromagnetic/,$(micromagnetic_objects))
