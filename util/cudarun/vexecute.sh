#!/bin/bash -i

# process commmand line arguments
while getopts d: option
do
case "${option}"
in
d) VDIR=${OPTARG};;
esac
done

# change into code directory
cd $VDIR

# Print working directory to screen
HN=`hostname`
LDIR=`pwd`
echo "Running code on $HN:$LDIR"

# get GPU info
nvidia-smi

# load shared libraries
module load gnu/6.3.0
module load cuda/11.0.3

# touch essential files to make sure they are updated via NFS
touch vampire-cuda
touch input
touch Co.mat

# execute vampire
nvprof --unified-memory-profiling off ./vampire-cuda
