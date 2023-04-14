#!/bin/bash

setting=${1}
cwd=`pwd`

echo  ${MOPSA_ROOT}/bin/evaluate3d $setting
time ${MOPSA_ROOT}/bin/evaluate3d $setting
