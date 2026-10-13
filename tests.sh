#!/bin/bash
#------------------------------------------------------
# simple bash script to build and run unit and
# integration tests
#------------------------------------------------------

usage(){
   echo "Usage: $0 [-a] [-i] [-u] [-h]"
   echo "   -a   run all tests (default if no option given)"
   echo "   -i   run integration tests"
   echo "   -u   run unit tests"
   echo "   -h   show this help"
}

integration=false
unit=false

while getopts "aiuh" flag
do
   case "${flag}" in
      a) integration=true; unit=true;;
      i) integration=true;;
      u) unit=true;;
      h) usage; exit 0;;
      *) usage; exit 1;;
   esac
done

# run everything if no test type is specified
if [ "$integration" = false ] && [ "$unit" = false ]; then
   integration=true
   unit=true
fi

# work relative to the location of this script
root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

status=0

if [ "$integration" = true ]; then
   if [ ! -x "$root/vampire-serial" ]; then
      echo "Error: $root/vampire-serial not found - build it first with 'make serial-llvm' (macOS) or 'make serial'"
      exit 1
   fi
   make -C "$root/test/integration" || exit 1
   cd "$root/test/integration" || exit 1
   time ./integration_tests || status=1
fi

if [ "$unit" = true ]; then
   make -C "$root/test/unit" || exit 1
   cd "$root/test/unit" || exit 1
   time ./unit_tests || status=1
fi

exit $status
