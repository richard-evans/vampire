#!/bin/bash

# terminal output colors
readonly green='\033[0;32m'
readonly red='\033[0;31m'
readonly nc='\033[0m'

function help_message {
    echo "Script to perform physical tests on VAMPIRE."
    echo "Author: S R H Morris"
    echo
    echo "Options:"
    echo " --help (-h): Prints this message."
    echo
    echo " --test (-t): Performs the test number supplied."
    echo "              If this option is not given all tests are performed."
    echo "              Valid numbers are 1 - Tests applied field and integrator."
    echo "                                2 - Tests anisotropy and thermal field."
    echo "                                3 - Tests exchange."
}

function cleanup {
    # restore old files
    mv output.bak output 2>/dev/null
    mv input.bak input   2>/dev/null
    mv Co.mat.bak Co.mat 2>/dev/null
}

function is_within_tolerance {
    local readonly error=$(echo $1 | sed 's/[eE]+*/\*10\^/')
    local readonly tolerance=$(echo $2 | sed 's/[eE]+*/\*10\^/')

    if (( $(echo $error '<' $tolerance | bc -l) )); then
        echo -e "${green}passed${nc} (maximum error $1)"
    else
        echo -e "${red}failed${nc} (maximum error $1)"
    fi
}

function applied_field {
    echo -n "Testing applied magnetic field..........."

    dir=tests/physical/AppliedField

    cp $dir/input input
    cp $dir/Co.mat Co.mat

    ./vampire &>/dev/null

    # check vampire output against analytic results
    $dir/applied_field_errors.py > applied_field_errors.dat
    max_error=$(grep "# maximum error = " applied_field_errors.dat | awk '{print $5}')
    is_within_tolerance $max_error 1e-6
}

function thermal {
    echo -n "Testing thermal effects.................."

    dir=tests/physical/Thermal

    cp $dir/input input
    cp $dir/Co.mat Co.mat

    ./vampire &>/dev/null
    max_error=$($dir/boltzmann_distribution.py)

    is_within_tolerance $max_error 0.01
}

function mag_vs_t {
    echo -n "Testing magnetisation with temperature..."

    dir=tests/physical/MvT

    cp $dir/input input
    cp $dir/Co.mat Co.mat

    ./vampire &>/dev/null
    max_error=$($dir/mvt.py)

    is_within_tolerance $max_error 0.01
}

function perform_test {

    case $1 in
        1)
            applied_field
            ;;
        2)
            thermal
            ;;
        3)
            mag_vs_t
            ;;
        *)
            echo -e "${red}Error: unknown test number $1. See --help for details."
            ;;
    esac
}


trap cleanup EXIT INT

# backup old files
mv output output.bak 2>/dev/null
mv input input.bak   2>/dev/null
mv Co.mat Co.mat.bak 2>/dev/null

while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
        -h|--help)
            help_message
            exit
            ;;
        -t|--test)
            perform_test $2
            shift 2
            ;;
        *)
            echo -e "${red}Error: unknown option $key. See --help for details."
            exit
            ;;
    esac
done
