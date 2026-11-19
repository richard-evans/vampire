#!/bin/bash
#------------------------------------------------------------------------------
#
#   This file is part of the VAMPIRE open source package under the
#   Free BSD licence (see licence file for details).
#
#------------------------------------------------------------------------------
# Chained HAMR simulation on the granular medium of hamr-granular-t0 with
# elliptical grain rounding, each stage started from the checkpoint of the
# previous one (sim:load-checkpoint=restart), so that the stages can run with
# different builds:
#
#   relax  Monte Carlo time series at 300 K, saving a checkpoint
#   write  zero temperature HAMR write restarted from the relaxed state, saving
#          a checkpoint and the final spins
#   load   zero step time series restarted from the end of the write, with the
#          hierarchical dipole field and its output (as in a read-back
#          simulation) and binary file-per-node configuration output: its spins
#          must be bit for bit those saved by the write
#
# Sourced by run_gpu_suite.sh and mgpu_check.sh (which set $results), from the
# repository root:
#
#   chain_run <stage> <label of previous stage or -> <label> <precision> <extra input line or ""> <command ...>
#
# Results go to $results/<label>/hamr-granular-t0/run, so that vtest.py compare
# and configs can be used on them.
#------------------------------------------------------------------------------
chain_case=developers/test/validation/cases/hamr-granular-t0

chain_run(){
	local stage=$1 src=$2 label=$3 pr=$4 extra=$5
	shift 5
	local d=$results/$label/hamr-granular-t0/run
	rm -rf $d; mkdir -p $d
	cp $chain_case/* $d/; rm $d/case.json
	[ $src != - ] && cp $results/$src/hamr-granular-t0/run/vampire*.chk $d/
	echo "{\"label\": \"$label\", \"precision\": \"$pr\"}" >$results/$label/meta.json
	{
		cat $chain_case/input
		echo
		echo "create:voronoi-elliptical-rounding = 0.5"
		case $stage in
			relax) printf '%s\n' "sim:program = time-series" "sim:integrator = monte-carlo" "sim:temperature = 300.0" \
				"sim:total-time-steps = 200" "sim:time-steps-increment = 1" "sim:save-checkpoint = end" ;;
			write) printf '%s\n' "sim:load-checkpoint = restart" "sim:save-checkpoint = end" "config:atoms = end" \
				"config:output-format = binary" ;;
			load) printf '%s\n' "sim:program = time-series" "sim:temperature = 0.0" "sim:total-time-steps = 0" \
				"sim:time-steps-increment = 1" "sim:load-checkpoint = restart" "config:atoms = end" "config:output-format = binary" \
				"config:output-mode = file-per-node" "config:output-nodes = 20" "dipole:solver = hierarchical" \
				"cells:macro-cell-size = 1.0 !nm" "cells:probe-non-magnetic-cells" "dipole:output-dipole-field" \
				"dipole:output-atomistic-dipole-field" ;;
		esac
		echo "$extra"
	} >$d/input
	(cd $d && "$@" >stdout 2>&1) || echo "chain $stage ($label): run failed (see $d/stdout)"
	[ $stage = relax ] || grep -q "Checkpoint file loaded" $d/log* || echo "chain $stage ($label): checkpoint not loaded (see $d/log)"
}
