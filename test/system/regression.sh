#!/bin/bash

RED='\033[0;31m'
GRN='\033[0;32m'
NCL='\033[0m' # No Color

# run tests
test_result=0
echo "==============================================================================="
echo " Running regression test suite for VAMPIRE"
echo "-------------------------------------------------------------------------------"
echo ""
#------------------------------------------------
# exchange tests
#------------------------------------------------
cd exchange_energy
bash test.sh
cd ..
echo $test_result

exit $test_result
