#!/bin/bash
while getopts u:d:n: option
do
case "${option}"
in
u) USER=${OPTARG};;
d) VDIR=${OPTARG};;
n) NODE=${OPTARG};;
esac
done

# log into compute node and execute code
ssh -t -q $USER@$NODE /bin/bash -l $VDIR/util/cudarun/vexecute.sh -d $VDIR
