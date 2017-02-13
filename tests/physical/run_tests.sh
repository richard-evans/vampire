#!/bin/bash

# backup input and material files
mv input input.bak
mv Co.mat Co.mat.bak

# use single spin in magnetic field input and material files
cp tests/physical/single_spin_input input
cp tests/physical/single_spin_mat Co.mat

./vampire &>/dev/null

# check vampire output against analytic results
tests/physical/single_spin_errors.py > single_spin_errors.dat

# restore old input and material file
mv input.bak input
mv Co.mat.bak Co.mat
