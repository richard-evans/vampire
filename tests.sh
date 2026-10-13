#!/bin/bash
#------------------------------------------------------
# simple bash script to build and run unit and
# integration tests
#------------------------------------------------------

usage(){
   echo "Usage: $0 [-a] [-i] [-e] [-u] [-h]"
   echo "   -a   run all tests, including extended integration tests"
   echo "   -i   run standard integration tests"
   echo "   -e   run extended (long running) integration tests"
   echo "   -u   run unit tests"
   echo "   -h   show this help"
   echo "With no option, the standard integration and unit tests are run."
}

integration=false
extended=false
unit=false

while getopts "aieuh" flag
do
   case "${flag}" in
      a) integration=true; extended=true; unit=true;;
      i) integration=true;;
      e) extended=true;;
      u) unit=true;;
      h) usage; exit 0;;
      *) usage; exit 1;;
   esac
done

# run standard integration and unit tests if no test type is specified
if [ "$integration" = false ] && [ "$extended" = false ] && [ "$unit" = false ]; then
   integration=true
   unit=true
fi

# work relative to the location of this script
root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

status=0

if [ "$integration" = true ] || [ "$extended" = true ]; then
   if [ ! -x "$root/vampire-serial" ]; then
      echo "Error: $root/vampire-serial not found - build it first with 'make serial-llvm' (macOS) or 'make serial'"
      exit 1
   fi
   make -C "$root/test/integration" || exit 1
   cd "$root/test/integration" || exit 1
   if [ "$integration" = true ] && [ "$extended" = true ]; then
      time ./integration_tests --all || status=1
   elif [ "$extended" = true ]; then
      time ./integration_tests --extended || status=1
   else
      time ./integration_tests || status=1
   fi
fi

if [ "$unit" = true ]; then
   make -C "$root/test/unit" || exit 1
   cd "$root/test/unit" || exit 1
   time ./unit_tests || status=1
fi

exit $status
