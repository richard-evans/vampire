#!/usr/bin/env python3
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
# Numeric comparison of two vampire output files for the regression tests.
#
# Usage: compare.py <expected> <obtained> [rel] [abs]
#
# Files are compared line by line and field by field. Numeric fields pass if
# |obtained - expected| <= abs + rel*|expected|, other fields must be equal.
# Lines starting with # are ignored. Exit code 0 on pass, 1 on fail.
#
import math
import sys

def data_lines(filename):
   with open(filename) as f:
      return [line.split() for line in f if line.strip() and not line.startswith("#")]

def to_float(s):
   try:
      return float(s)
   except ValueError:
      return None

def main():

   if len(sys.argv) not in (3, 5):
      print("Usage: compare.py <expected> <obtained> [rel] [abs]")
      return 2

   expected = data_lines(sys.argv[1])
   obtained = data_lines(sys.argv[2])
   rel = float(sys.argv[3]) if len(sys.argv) == 5 else 1.0e-12
   tol = float(sys.argv[4]) if len(sys.argv) == 5 else 1.0e-14

   if len(expected) != len(obtained):
      print("different number of lines: expected {}, obtained {}".format(len(expected), len(obtained)))
      return 1

   # report the worst field and the number of failing fields
   num_fail = 0
   worst = None
   for line, (e_fields, o_fields) in enumerate(zip(expected, obtained), start=1):
      if len(e_fields) != len(o_fields):
         print("line {}: different number of fields: expected {}, obtained {}".format(line, len(e_fields), len(o_fields)))
         return 1
      for col, (e, o) in enumerate(zip(e_fields, o_fields), start=1):
         ev = to_float(e)
         ov = to_float(o)
         if ev is None or ov is None:
            if e != o:
               num_fail += 1
               if worst is None: worst = (float("inf"), line, col, e, o)
            continue
         # nan and inf only match themselves
         if not (math.isfinite(ev) and math.isfinite(ov)):
            if not (ev == ov or (math.isnan(ev) and math.isnan(ov))):
               num_fail += 1
               if worst is None or worst[0] != float("inf"): worst = (float("inf"), line, col, e, o)
            continue
         diff = abs(ov - ev)
         bound = tol + rel*abs(ev)
         if diff > bound:
            excess = diff/bound if bound > 0.0 else float("inf")
            num_fail += 1
            if worst is None or excess > worst[0]: worst = (excess, line, col, e, o)

   if num_fail > 0:
      excess, line, col, e, o = worst
      print("{} fields differ; worst at line {} column {}: expected {}, obtained {} ({:.3g} x tolerance)".format(
            num_fail, line, col, e, o, excess))
      return 1

   return 0

if __name__ == "__main__":
   sys.exit(main())
