#!/bin/bash
#------------------------------------------------------------------------------
#
#   This file is part of the VAMPIRE open source package under the
#   Free BSD licence (see licence file for details).
#
#------------------------------------------------------------------------------
# Incremental builds of the working tree for the test scripts (sourced, not
# run). Each build has its own directory, a copy of makefile, src, hdr and obj,
# so that builds of different targets and precisions keep their objects and the
# working tree is not cleaned.
#
#    incremental_build <directory> <target> [make variables ...]
#
# recompiles only the source files that changed since the last build in that
# directory, or everything if a header, a makefile or the make variables
# changed: the makefile does not track header dependencies. Run from the
# repository root; JOBS sets the number of parallel jobs.
#------------------------------------------------------------------------------
command -v rsync > /dev/null || { echo "rsync is needed for the incremental builds"; exit 1; }

incremental_build(){
	local b=$1 target=$2
	shift 2
	local stamp="$target $* ${MAKEFLAGS:-}" f stem
	mkdir -p $b
	local changed=$(rsync -a --delete --exclude='*.o' --out-format='%n' makefile src hdr obj $b/)
	# objects of changed sources are removed, so make rebuilds them whatever the file times
	for f in $(echo "$changed" | grep -E '^src/.*\.(cpp|cu)$'); do
		stem=${f#src/}; stem=${stem%.*}
		rm -f $b/obj/$stem.o $b/obj/${stem}_*.o
	done
	if [ "$(cat $b/build-settings.txt 2> /dev/null)" != "$stamp" ] || echo "$changed" | grep -qE '\.(h|hpp|cuh|inc)$|makefile$'; then
		make -C $b clean > /dev/null
	fi
	echo "$stamp" > $b/build-settings.txt
	make -C $b -j${JOBS:-$(nproc)} $target "$@"
}
