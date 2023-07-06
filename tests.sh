#!/bin/bash
#------------------------------------------------------
# simple bash script to run unit and integration tests
#------------------------------------------------------

while getopts ::aiu flag
do
    case "${flag}" in
        a) all=true;;
        i) integration=true;;
        u) unit=true;;
    esac
done

if [ "$all" = true ]; then
    integration=true
    unit=true
fi

if [ "$integration" = true ]; then
    echo "====================================================================="
    echo "      Running integration tests"
    echo "====================================================================="
    cd test/integration/
    time ./integration_tests
    cd ../../
fi
if [ "$unit" = true ]; then
    echo "Running unit tests"
fi
