#!/bin/bash
#------------------------------------------------------------------------------
#
#   This file is part of the VAMPIRE open source package under the
#   Free BSD licence (see licence file for details).
#
#------------------------------------------------------------------------------
# Checks the macro-cell dipole solvers on a uniformly magnetised Co film
# (40 x 12 x 3 nm, 1 nm cells, longer than the range over which the
# hierarchical solver sums cell by cell):
#
#   - dipole:on-the-fly-tensors gives fields identical to stored tensors
#   - the hierarchical field agrees with the tensor field (largest difference
#     below TOL tesla, default 1e-3)
#
#    developers/test/validation/dipole_check.sh [vampire binary] [MPI ranks]
#
# The binary defaults to ./vampire-serial; with a rank count the runs use
# mpirun -np <ranks>. Results go to a temporary directory, which is removed
# unless KEEP=1.
#------------------------------------------------------------------------------
set -eu

exe=$(realpath "${1:-./vampire-serial}")
ranks=${2:-}
tol=${TOL:-1e-3}
work=$(mktemp -d)
[ "${KEEP:-0}" = 1 ] && echo "results in $work" || trap 'rm -rf "$work"' EXIT

cat > "$work/Co.mat" << 'EOF'
material:num-materials=1
material[1]:material-name=Co
material[1]:damping-constant=1.0
material[1]:exchange-matrix[1]=11.2e-21
material[1]:atomic-spin-moment=1.72 !muB
material[1]:material-element=Ag
material[1]:initial-spin-direction=0,0,1
EOF

run(){ # name solver [extra input line]
   mkdir "$work/$1"
   cp "$work/Co.mat" "$work/$1/"
   cat > "$work/$1/input" << EOF
dimensions:unit-cell-size = 3.54 !A
dimensions:system-size-x = 40 !nm
dimensions:system-size-y = 12 !nm
dimensions:system-size-z = 3 !nm
material:file = Co.mat
cells:macro-cell-size = 1 !nm
dipole:solver = $2
dipole:output-dipole-field
${3:-}
sim:temperature = 0.0
sim:time-step = 1 !fs
sim:total-time-steps = 0
sim:program = time-series
sim:integrator = llg-heun
output:real-time
EOF
   if [ -n "$ranks" ]; then
      (cd "$work/$1" && mpirun -np "$ranks" "$exe" > stdout 2>&1) || { tail -n 20 "$work/$1/stdout"; exit 1; }
   else
      (cd "$work/$1" && "$exe" > stdout 2>&1) || { tail -n 20 "$work/$1/stdout"; exit 1; }
   fi
}

run tensor tensor
run tensor-otf tensor dipole:on-the-fly-tensors
run hierarchical hierarchical

fail=0
if cmp -s "$work/tensor/dipole-field.txt" "$work/tensor-otf/dipole-field.txt"; then
   echo "PASS tensor: on-the-fly tensors give identical fields"
else
   echo "FAIL tensor: on-the-fly tensors give different fields"; fail=1
fi

python3 - "$work" "$tol" << 'EOF' || fail=1
import sys
import numpy as np
work, tol = sys.argv[1], float(sys.argv[2])
t = np.loadtxt(work + "/tensor/dipole-field.txt")
h = np.loadtxt(work + "/hierarchical/dipole-field.txt")
diff = np.abs(h[:, 6:9] - t[:, 6:9]).max()
ok = diff < tol
print(f"{'PASS' if ok else 'FAIL'} hierarchical vs tensor: largest field difference {diff:.2e} T (tolerance {tol:g} T, "
      f"largest field {np.abs(t[:, 6:9]).max():.3f} T)")
sys.exit(0 if ok else 1)
EOF

exit $fail
