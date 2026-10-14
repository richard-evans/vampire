#!/bin/bash
#---------------------------------------------------------------
# A simple script to get (most) of the vampire input parameters
#---------------------------------------------------------------
grep -E 'test=|test =' src/*/interface.cpp src/vio/match.cpp | tr -d '\"; '  | tr '=' ' ' | awk '{print $2}' > input_paramaters.txt
grep -E 'test=|test =' src/*/interface.cpp src/vio/match.cpp | tr -d '\"; '  | tr '=/' ' ' | awk '{print $2}' > module.txt

grep -A7 -E 'test=|test =' src/*/interface.cpp src/vio/match.cpp | grep -E "test=|test =|double|int|bool|check|string " > extra_info.txt


