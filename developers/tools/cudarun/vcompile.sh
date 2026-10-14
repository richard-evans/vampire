#!/bin/bash -i
# ^ interactive shell needed for module command

while getopts d: option
do
case "${option}"
in
d) VDIR=${OPTARG};;
esac
done

# change into source code directory
cd $VDIR

# Output debugging information
#pwd
#hostname

# load shared libraries
module load gnu/6.3.0
module load cuda/11.0.3

# compile code
make cuda -j 16


