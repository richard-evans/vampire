#!/bin/bash
#------------------------------------------------------
# simple bash script to build and run unit,
# integration and regression tests
#------------------------------------------------------

usage(){
   echo "Usage: $0 [-a] [-i] [-e] [-u] [-r] [-v] [-m [N]] [-b] [-c] [-s] [-h]"
   echo "   -a   run all tests, including extended integration tests, regression tests and the style check"
   echo "   -i   run standard integration tests"
   echo "   -e   run extended (long running) integration tests"
   echo "   -u   run unit tests"
   echo "   -r   run regression tests with vampire-serial"
   echo "   -v   run the physics validation suite with vampire-serial against the CPU references"
   echo "   -m   run regression tests with vampire-parallel on N processes (default 2)"
   echo "   -b   build matrix: clean build of every build variant in a temporary copy"
   echo "   -c   coverage report from a debug build in a temporary copy (needs gcovr)"
   echo "   -s   check source code style (tabs, trailing whitespace, licence headers)"
   echo "   -h   show this help"
   echo "With no option, the standard integration and unit tests are run."
}

integration=false
extended=false
unit=false
regression=false
validation=false
mpi_processes=0
build=false
coverage=false
style=false

while getopts "aieurvmbcsh" flag
do
   case "${flag}" in
      a) integration=true; extended=true; unit=true; regression=true; validation=true; style=true;;
      i) integration=true;;
      e) extended=true;;
      u) unit=true;;
      r) regression=true;;
      v) validation=true;;
      m) mpi_processes=2
         # optional number of processes as the next argument
         if [[ "${!OPTIND}" =~ ^[1-9][0-9]*$ ]]; then
            mpi_processes="${!OPTIND}"
            OPTIND=$((OPTIND+1))
         fi;;
      b) build=true;;
      c) coverage=true;;
      s) style=true;;
      h) usage; exit 0;;
      *) usage; exit 1;;
   esac
done

# run standard integration and unit tests if no test type is specified
if [ "$integration" = false ] && [ "$extended" = false ] && [ "$unit" = false ] && [ "$regression" = false ] && [ "$validation" = false ] &&
   [ "$mpi_processes" -eq 0 ] && [ "$build" = false ] && [ "$coverage" = false ] && [ "$style" = false ]; then
   integration=true
   unit=true
fi

# work relative to the vampire root, the parent of this script's directory
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

status=0

if [ "$style" = true ]; then
   # print errors and the summary; warnings are listed by running the script directly
   style_output="$("$root/developers/tools/check_style.sh")" || status=1
   echo "$style_output" | grep -v ": warning: "
fi

#------------------------------------------------------
# copy the working tree (tracked and untracked files,
# without ignored files such as objects) to a new
# temporary directory and print its name
#------------------------------------------------------
copy_tree(){
   local dir
   dir="$(mktemp -d "${TMPDIR:-/tmp}/vampire-$1.XXXXXX")"
   ( cd "$root" && git ls-files -co --exclude-standard -z | xargs -0 tar cf - ) | tar xf - -C "$dir"
   echo "$dir"
}

if [ "$integration" = true ] || [ "$extended" = true ]; then
   if [ ! -x "$root/vampire-serial" ]; then
      echo "Error: $root/vampire-serial not found - build it first with 'make serial-llvm' (macOS) or 'make serial'"
      exit 1
   fi
   make -C "$root/developers/test/integration" || exit 1
   cd "$root/developers/test/integration" || exit 1
   if [ "$integration" = true ] && [ "$extended" = true ]; then
      time ./integration_tests --all || status=1
   elif [ "$extended" = true ]; then
      time ./integration_tests --extended || status=1
   else
      time ./integration_tests || status=1
   fi
fi

if [ "$unit" = true ]; then
   make -C "$root/developers/test/unit" || exit 1
   cd "$root/developers/test/unit" || exit 1
   time ./unit_tests || status=1
fi

if [ "$regression" = true ]; then
   time "$root/developers/test/regression/run_regression.sh" --serial || status=1
fi

if [ "$validation" = true ]; then
   cd "$root" || exit 1
   python3 developers/test/validation/vtest.py run --exe ./vampire-serial --label serial --jobs 4 || status=1
   python3 developers/test/validation/vtest.py compare --test serial --ref cpu || status=1
fi

if [ "$mpi_processes" -gt 0 ]; then
   time "$root/developers/test/regression/run_regression.sh" --mpi "$mpi_processes" || status=1
fi

#------------------------------------------------------
# Build matrix
#------------------------------------------------------
if [ "$build" = true ]; then

   # build variants as "name:make arguments"
   if [ "$(uname -s)" = Darwin ]; then
      variants=("serial-llvm:serial-llvm" "parallel-llvm:parallel-llvm")
   else
      variants=("serial:serial" "parallel:parallel")
   fi
   variants+=("serial-debug:serial-debug" "parallel-debug:parallel-debug" "vdc:vdc")

   # FFT build if FFTW is installed
   for prefix in /opt/local /usr/local /opt/homebrew /usr; do
      if [ -f "$prefix/include/fftw3.h" ]; then
         variants+=("${variants[0]%%:*}-fft:${variants[0]#*:} FFTW=-DFFT\ -I$prefix/include LIBS=-lm\ -lfftw3\ -L$prefix/lib")
         break
      fi
   done

   echo "---------------------------------------------------------------------"
   echo "    Build matrix"
   echo "---------------------------------------------------------------------"
   for variant in "${variants[@]}"; do
      name="${variant%%:*}"
      dir="$(copy_tree build-$name)"
      printf "%-30s : " "$name"
      start=$SECONDS
      eval "make -C \"$dir\" -j 8 ${variant#*:}" > "$dir/build.log" 2>&1
      result=$?
      warnings=$(grep -c "warning:" "$dir/build.log")
      if [ $result -eq 0 ]; then
         echo "OK ($warnings warnings, $((SECONDS-start)) s)"
         rm -rf "$dir"
      else
         echo "FAIL (see $dir/build.log)"
         status=1
      fi
   done
   echo "---------------------------------------------------------------------"

fi

#------------------------------------------------------
# Coverage
#------------------------------------------------------
if [ "$coverage" = true ]; then

   # gcovr may be installed as a python module only (pip install --user gcovr)
   if command -v gcovr > /dev/null; then
      gcovr="gcovr"
   elif python3 -m gcovr --version > /dev/null 2>&1; then
      gcovr="python3 -m gcovr"
   else
      echo "Error: coverage needs gcovr (for example 'sudo port install py-gcovr' or 'pip install gcovr')"
      exit 1
   fi

   # debug build instrumented for coverage, used in place of vampire-serial by the tests
   dir="$(copy_tree coverage)"
   echo "Building coverage executable in $dir"
   make -C "$dir" -j 8 serial-debug > "$dir/build.log" 2>&1 || { echo "Error: build failed, see $dir/build.log"; exit 1; }
   ln -s vampire-serial-debug "$dir/vampire-serial"

   make -C "$dir/developers/test/integration" > /dev/null || exit 1
   ( cd "$dir/developers/test/integration" && ./integration_tests --all )
   "$dir/developers/test/regression/run_regression.sh" --serial

   mkdir -p "$dir/coverage"
   ( cd "$dir" && $gcovr --root . --filter src/ --exclude src/qvoronoi/ --object-directory obj \
        --html-details coverage/index.html --txt coverage/summary.txt ) || exit 1
   cat "$dir/coverage/summary.txt"
   echo "Coverage report: $dir/coverage/index.html"

fi

exit $status
