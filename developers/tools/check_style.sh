#!/bin/bash
#------------------------------------------------------------------------------
#
#   This file is part of the VAMPIRE open source package under the
#   Free BSD licence (see licence file for details).
#
#   (c) R F L Evans 2026. All rights reserved.
#
#------------------------------------------------------------------------------
#
#   Checks the C++, CUDA and OpenCL sources in src/ and hdr/ for the vampire
#   coding style:
#
#      - hard tabs (indentation is 3 spaces)
#      - trailing whitespace
#      - CRLF line endings
#      - missing newline at the end of the file
#      - missing licence header (the standard block produced by
#        developers/tools/initialise_new_vampire_module.cpp, Free BSD or GNU GPL v2)
#      - includes of another module's internal.hpp ("../<module>/internal.hpp")
#
#   Cross-module internal.hpp includes are reported as warnings unless
#   --strict is given, until the remaining ones have been removed.
#
#   Usage: developers/tools/check_style.sh [--strict] [file...]
#
#   With no files, all tracked sources are checked. The exit status is 1 if
#   any error is found.
#
#------------------------------------------------------------------------------

usage(){
   echo "Usage: $0 [--strict] [file...]"
   echo "   --strict   treat cross-module internal.hpp includes as errors"
   echo "   -h         show this help"
}

strict=false
while [ $# -gt 0 ]; do
   case "$1" in
      --strict) strict=true; shift;;
      -h|--help) usage; exit 0;;
      -*) usage; exit 1;;
      *) break;;
   esac
done

# work from the repository root
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$root" || exit 1

# vendored third-party files that keep their own formatting and licence
vendored='^(hdr/mtrand\.hpp|src/random/mtrand\.cpp|src/cuda/cuda_timer\.h)$'

if [ $# -gt 0 ]; then
   files=("$@")
else
   files=()
   while IFS= read -r f; do
      files+=("$f")
   done < <(git ls-files 'src/*.cpp' 'src/*.hpp' 'src/*.h' 'src/*.cu' 'src/*.cuh' 'src/*.cl' 'hdr/*.hpp' 'hdr/*.h' |
            grep -v -E "$vendored")
fi

if [ ${#files[@]} -eq 0 ]; then
   echo "No files to check"
   exit 0
fi

# line-by-line checks; prints "file:line: message" and a summary count per check
awk -v strict="$strict" '
   function header_done(){
      if(file != "" && !licence) report(file, 1, "licence", "missing licence header")
   }
   function report(f, n, kind, msg){
      count[kind]++
      if(kind == "include" && strict != "true"){
         printf "%s:%d: warning: %s\n", f, n, msg
      }
      else {
         printf "%s:%d: error: %s\n", f, n, msg
         errors++
      }
   }
   FNR == 1 { header_done(); file = FILENAME; licence = 0; nfiles++ }
   FNR <= 30 && /^\/\/   (Free BSD|GNU GPL \(version 2\)) licence \(see licence file for details\)\.$/ { licence = 1 }
   /\r$/                 { report(FILENAME, FNR, "crlf", "CRLF line ending") }
   /\t/                  { report(FILENAME, FNR, "tab", "hard tab") }
   /[ \t]+\r?$/          { report(FILENAME, FNR, "trailing", "trailing whitespace") }
   /^[ \t]*#[ \t]*include[ \t]*"\.\.\/[^"]*internal\.hpp"/ {
      report(FILENAME, FNR, "include", "includes another module'"'"'s internal.hpp")
   }
   END {
      header_done()
      printf "Checked %d files: %d errors", nfiles, errors
      printf " (tabs %d, trailing whitespace %d, CRLF %d, licence %d, cross-module includes %d)\n",
             count["tab"], count["trailing"], count["crlf"], count["licence"], count["include"]
      exit (errors > 0)
   }
' "${files[@]}"
status=$?

# missing final newline (awk cannot see this)
for f in "${files[@]}"; do
   if [ -s "$f" ] && [ -n "$(tail -c 1 "$f")" ]; then
      echo "$f: error: no newline at end of file"
      status=1
   fi
done

exit $status
