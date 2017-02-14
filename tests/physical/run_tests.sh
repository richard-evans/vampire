#!/bin/bash

# terminal output colors
readonly green='\033[0;32m'
readonly red='\033[0;31m'
readonly nc='\033[0m'


function is_within_tolerance {
    local readonly tolerance=0.000001
    local readonly error=$(echo $1 | sed 's/[eE]+*/\*10\^/')

    if (( $(echo $error '<' $tolerance | bc -l) )); then
	echo -e "${green}passed${nc} (maximum error $1)"
    else
	echo -e "${red}failed${nc} (maximum error $1)"
    fi
}


# backup input and material files
mv input input.bak
mv Co.mat Co.mat.bak



#### Test One ####
echo -n "Testing applied magnetic field....."

dir=tests/physical/AppliedField

# use single spin in magnetic field input and material files
cp $dir/input input
cp $dir/Co.mat Co.mat

./vampire &>/dev/null

# check vampire output against analytic results
$dir/applied_field_errors.py > applied_field_errors.dat
max_error=$(grep "# maximum error = " applied_field_errors.dat | awk '{print $5}')
is_within_tolerance $max_error
#### End Test One ####

#### Test Two ####
echo -n "Testing thermal effects............"

dir=tests/physical/Thermal

cp $dir/input input
cp $dir/Co.mat Co.mat

./vampire

is_within_tolerance 0.01


# restore old input and material file
mv input.bak input
mv Co.mat.bak Co.mat
