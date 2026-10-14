#!/bin/bash
#------------------------------------------------------------------------------
#
#   This file is part of the VAMPIRE open source package under the
#   Free BSD licence (see licence file for details).
#
#   (c) Richard F L Evans 2026. All rights reserved.
#
#   Email: richard.evans@york.ac.uk
#
#------------------------------------------------------------------------------
#
# Golden regression tests: runs each case in cases/<module>/<case>/ and
# compares the output file against the stored golden output.
#
# A case directory contains:
#    input, *.mat, [*.ucf, ...]   inputs, copied into a temporary directory
#    expected-serial.dat          golden output (comment lines stripped)
#    expected-mpi-np<N>.dat       golden output for vampire-parallel on N processes
#    tolerance                    optional, "rel <r> abs <a>" (default rel 1e-12 abs 1e-14)
#    files                        optional, other output files to compare, one per line (glob patterns
#                                 allowed); their golden files are expected-serial.<file> and
#                                 expected-mpi-np<N>.<file>. Binary files are compared byte for byte.
#    expected-error               for input that must be rejected: vampire must exit with an error and
#                                 its terminal output must contain each line of this file (all modes)
#
# A case is only run in a mode for which it has a golden file, unless --bless is given.
# Each run is stopped after TIME_LIMIT seconds (default 300), so a hang is reported as a failure.
#
#------------------------------------------------------------------------------

usage(){
   echo "Usage: $0 [--serial | --mpi N] [--exact] [--bless] [--case <glob>] [--keep]"
   echo "   --serial       run vampire-serial (default)"
   echo "   --mpi N        run vampire-parallel on N processes"
   echo "   --exact        compare byte for byte instead of with a numerical tolerance"
   echo "   --bless        write the obtained output as the new golden file"
   echo "   --case <glob>  only run cases matching <module>/<case> glob, eg 'integrate/llg-*'"
   echo "   --keep         keep the working directory of every case, not just failed ones"
}

mode=serial
nprocs=""
exact=false
bless=false
pattern="*/*"
keep=false

while [ $# -gt 0 ]; do
   case "$1" in
      --serial) mode=serial;;
      --mpi)    mode=mpi; nprocs="$2"; shift;;
      --exact)  exact=true;;
      --bless)  bless=true;;
      --case)   pattern="$2"; shift;;
      --keep)   keep=true;;
      -h|--help) usage; exit 0;;
      *) usage; exit 1;;
   esac
   shift
done

if [ "$mode" = mpi ] && ! [[ "$nprocs" =~ ^[1-9][0-9]*$ ]]; then
   echo "Error: --mpi needs a number of processes"
   exit 1
fi

# work relative to the location of this script
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../../.." && pwd)"

if [ "$mode" = serial ]; then
   executable="$root/vampire-serial"
   golden=expected-serial.dat
   run="$executable"
else
   executable="$root/vampire-parallel"
   golden=expected-mpi-np$nprocs.dat
   run="${MPIRUN:-mpirun} -np $nprocs $executable"
fi

if [ ! -x "$executable" ]; then
   echo "Error: $executable not found - build it first"
   exit 1
fi

time_limit="${TIME_LIMIT:-300}"

work="$(mktemp -d "${TMPDIR:-/tmp}/vampire-regression.XXXXXX")"

pass=0
fail=0
skip=0
blessed=0

echo "---------------------------------------------------------------------"
echo "    Running regression tests ($mode${nprocs:+, np = $nprocs})"
echo "---------------------------------------------------------------------"

cd "$here/cases" || exit 1
for name in $pattern; do

   [ -f "$name/input" ] || continue

   # skip cases without a golden file for this mode
   if [ "$bless" = false ] && [ ! -f "$name/$golden" ] && [ ! -f "$name/expected-error" ]; then
      skip=$((skip+1))
      continue
   fi

   printf "%-60s : " "$name"

   # run in a clean copy of the case directory
   dir="$work/$name"
   mkdir -p "$dir"
   for f in "$name"/*; do
      case "$(basename "$f")" in
         expected-*|tolerance|files) ;;
         *) cp -R "$f" "$dir/";;
      esac
   done

   start=$SECONDS
   ( cd "$dir" && perl -e 'alarm shift; exec @ARGV' "$time_limit" $run > stdout 2>&1 ) 2> /dev/null
   status=$?
   elapsed=$((SECONDS-start))

   if [ $status -eq 142 ]; then
      echo "FAIL (timed out after $time_limit s, see $dir/stdout)"
      fail=$((fail+1))
      continue
   fi

   # input that must be rejected with an error message
   if [ -f "$name/expected-error" ]; then
      missing=""
      while IFS= read -r text; do
         [ -n "$text" ] && ! grep -qF -- "$text" "$dir/stdout" && missing="$missing \"$text\""
      done < "$name/expected-error"
      if [ $status -eq 0 ]; then
         echo "FAIL (vampire accepted the input, see $dir/stdout)"
         fail=$((fail+1))
      elif [ -n "$missing" ]; then
         echo "FAIL (error message missing$missing, see $dir/stdout)"
         fail=$((fail+1))
      else
         echo "OK (rejected, ${elapsed} s)"
         pass=$((pass+1))
         [ "$keep" = true ] || rm -rf "$dir"
      fi
      continue
   fi

   if [ $status -ne 0 ]; then
      echo "FAIL (vampire exited with status $status, see $dir/stdout)"
      fail=$((fail+1))
      continue
   fi

   # some programs never write the output file, so treat a missing one as empty
   [ -f "$dir/output" ] || : > "$dir/output"
   # output files to compare, as pairs of obtained file and golden file. Entries in the files list may be
   # glob patterns, since the files written can differ between serial and mpi runs.
   result=""
   prefix="$name/${golden%.dat}"
   compare_list=("output" "$name/$golden")
   if [ -f "$name/files" ]; then
      while read -r pattern_f; do
         [ -n "$pattern_f" ] || continue
         matches="$(cd "$dir" && ls -1d $pattern_f 2> /dev/null)"
         [ -n "$matches" ] || result="$result $pattern_f: not written;"
         for f in $matches; do
            compare_list+=("$f" "$prefix.$f")
         done
      done < "$name/files"
   fi

   if [ "$bless" = true ]; then
      # remove goldens of this mode from earlier blessings before writing new ones
      for g in "$prefix".*; do
         [ "$g" = "$name/$golden" ] || rm -f "$g"
      done
   else
      # every golden file of this mode must have been written again
      for g in "$prefix".*; do
         [ -f "$g" ] && [ "$g" != "$name/$golden" ] || continue
         f="${g#$prefix.}"
         [ -f "$dir/$f" ] || result="$result $f: not written;"
      done
   fi

   for ((i=0; i<${#compare_list[@]}; i+=2)); do
      obtained="${compare_list[i]}"
      expected="${compare_list[i+1]}"
      if [ ! -f "$dir/$obtained" ]; then
         result="$result $obtained: not written;"
         continue
      fi
      # binary files are compared byte for byte, text files without comment lines
      binary=false
      if [ -s "$dir/$obtained" ] && ! grep -qI . "$dir/$obtained"; then
         binary=true
         cp "$dir/$obtained" "$dir/$obtained.obtained"
      else
         grep -v "^#" "$dir/$obtained" > "$dir/$obtained.obtained"
      fi
      if [ "$bless" = true ]; then
         cp "$dir/$obtained.obtained" "$expected"
      elif [ ! -f "$expected" ]; then
         result="$result $obtained: no golden file;"
      elif [ "$exact" = true ] || [ "$binary" = true ]; then
         cmp -s "$expected" "$dir/$obtained.obtained" || result="$result $obtained: files differ;"
      else
         tolerance=""
         [ -f "$name/tolerance" ] && tolerance="$(awk '{for(i=1;i<NF;i++){if($i=="rel")r=$(i+1);if($i=="abs")a=$(i+1)}} END{print r, a}' "$name/tolerance")"
         diffs="$(python3 -I "$here/compare.py" "$expected" "$dir/$obtained.obtained" $tolerance)" || result="$result $obtained: ${diffs:-compare failed};"
      fi
   done

   if [ -z "$result" ] && [ "$bless" = true ]; then
      echo "BLESSED ($(( ${#compare_list[@]}/2 )) files, ${elapsed} s)"
      blessed=$((blessed+1))
      [ "$keep" = true ] || rm -rf "$dir"
   elif [ -z "$result" ]; then
      echo "OK (${elapsed} s)"
      pass=$((pass+1))
      [ "$keep" = true ] || rm -rf "$dir"
   else
      echo "FAIL |$result see $dir"
      fail=$((fail+1))
   fi

done

# remove the working directory unless it holds failed or kept cases
[ -z "$(find "$work" -type f)" ] && rm -rf "$work"

echo "---------------------------------------------------------------------"
if [ "$bless" = true ]; then
   echo "Blessed $blessed cases, $fail failed to run"
else
   echo "Passed $pass, failed $fail, skipped $skip (no $golden)"
   if [ $fail -gt 0 ]; then echo "OVERALL FAIL"; else echo "OVERALL PASS"; fi
fi
echo "---------------------------------------------------------------------"

[ $fail -eq 0 ]
