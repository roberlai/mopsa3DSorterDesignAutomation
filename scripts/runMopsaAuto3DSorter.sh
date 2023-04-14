#!/bin/bash

comsol server -silent -port 2036 -login never &
sleep 10
matlab_setting=${1}

dir=$(dirname $(realpath $matlab_setting))
filename=$(basename -s .m $matlab_setting)

cd $MOPSA_ROOT/mopsa_auto/3dsorter
echo matlab -nodesktop -r "addpath('$dir'); rehash toolboxcache; declareGlobalVar, $filename, main_from_given_recipe, quit"
time matlab -nodesktop -r "addpath('$dir'); rehash toolboxcache; declareGlobalVar, $filename, main_from_given_recipe, quit"
