#!/bin/bash

setting=${1}
cwd=`pwd`

cd ${MOPSA_ROOT}/build
make -j16

reval=$?
if [ $reval -ne 0 ]; then
  echo "Fail to build"
  exit
fi

cd $cwd

echo ${MOPSA_ROOT}/bin/main3d $setting
time ${MOPSA_ROOT}/bin/main3d $setting
