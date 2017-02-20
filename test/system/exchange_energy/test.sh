#!/bin/bash

#vampire executable
vampire='../../../../vampire'

RED='\033[0;31m'
GRN='\033[0;32m'
NCL='\033[0m' # No Color

# run tests
test_result=0
echo "-------------------------------------------------------------------------------"
echo " Exchange tests"
echo "-------------------------------------------------------------------------------"
#------------------------------------------------
# simple cubic
#------------------------------------------------
cd sc
$vampire $1 > /dev/null
tail -n1 output > result
difference=`diff result expected`
if [ "$difference" != "" ]
then
   #       |----------------------------------------|---------|
   echo -e "Simple cubic:                           [${RED}failed${NCL}]"
   echo "   Expected | Result"
   paste expected result
   test_result=1
else
   #       |----------------------------------------|---------|
   echo -e "Simple cubic:                           [${GRN}OK${NCL}]"
fi
cd ..

#------------------------------------------------
# body centred cubic
#------------------------------------------------
cd bcc
$vampire $1 > /dev/null
tail -n1 output > result
difference=`diff result expected`
if [ "$difference" != "" ]
then
   #       |----------------------------------------|---------|
   echo -e "Body centred cubic:                     [${RED}failed${NCL}]"
   echo "   Expected | Result"
   paste expected result
   test_result=1
else
   #       |----------------------------------------|---------|
   echo -e "Body centred cubic:                     [${GRN}OK${NCL}]"
fi
cd ..

#------------------------------------------------
# face centred cubic
#------------------------------------------------
cd fcc
$vampire $1 > /dev/null
tail -n1 output > result
difference=`diff result expected`
if [ "$difference" != "" ]
then
   #       |----------------------------------------|---------|
   echo -e "Face centred cubic:                     [${RED}failed${NCL}]"
   echo "   Expected | Result"
   paste expected result
   test_result=1
else
   #       |----------------------------------------|---------|
   echo -e "Face centred cubic:                     [${GRN}OK${NCL}]"
fi
cd ..

#------------------------------------------------
# hexagonal close packed
#------------------------------------------------
cd hcp
$vampire $1 > /dev/null
tail -n1 output > result
difference=`diff result expected`
if [ "$difference" != "" ]
then
   #       |----------------------------------------|---------|
   echo -e "Hexagonal close packed:                 [${RED}failed${NCL}]"
   echo "   Expected | Result"
   paste expected result
   test_result=1
else
   #       |----------------------------------------|---------|
   echo -e "Hexagonal close packed:                 [${GRN}OK${NCL}]"
fi
cd ..

exit $test_result
