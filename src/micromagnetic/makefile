#--------------------------------------------------------------
#          Makefile for micromagnetic module
#--------------------------------------------------------------

# List module object filenames
micromagnetic_objects =\
atomistic_llg_heun.o \
bias_magnets.o \
calculate_A.o \
calculate_alpha.o \
calculate_chi_para.o \
calculate_chi_perp.o \
calculate_gamma.o \
calculate_ku.o \
calculate_ms.o \
calculate_tc.o \
calculate_stt.o \
data.o \
initialize.o \
interface.o \
fieldsllb.o \
fieldsllg.o \
micromagnetic_llb.o \
micromagnetic_llg.o \
print.o \
resistance.o \
output.o
# Append module objects to global tree
OBJECTS+=$(addprefix obj/micromagnetic/,$(micromagnetic_objects))
