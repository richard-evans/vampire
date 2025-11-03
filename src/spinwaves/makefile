#--------------------------------------------------------------
#          Makefile for sw module
#--------------------------------------------------------------

# List module object filenames
spinwave_objects =\
data.o \
built_in_paths.o \
paths.o \
initialize.o \
interface.o \
fft_in_time.o \
fft_in_time_options.o \
fft_in_space.o \
util.o

# Append module objects to global tree
OBJECTS+=$(addprefix obj/spinwaves/,$(spinwave_objects))
