#!/bin/bash
#Do not forget to install libgfortran3  on your system!
module load releases/2020b
module load CMake/3.18.4-GCCcore-10.2.0 
module load Lua/5.4.2-GCCcore-10.2.0
module load GCC/10.2.0
module load CGAL/5.2-gompi-2020b 

cd ../../../

if [ ! -d "dependencies" ]; then
    mkdir dependencies
fi

cd dependencies/

if [ ! -d "gmsh-4.8.3-Linux64-sdk" ]; then
  wget http://gmsh.info/bin/Linux/gmsh-4.8.3-Linux64-sdk.tgz
  tar -xf gmsh-4.8.3-Linux64-sdk.tgz 
  rm -rf gmsh-4.8.3-Linux64-sdk.tgz 
fi

if [ ! -d "eigen-3.3.9" ]; then
  wget https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.gz
  tar -xf eigen-3.3.9.tar.gz
  rm -rf eigen-3.3.9.tar.gz
fi

if [ ! -d "sol" ]; then
  mkdir sol
  cd sol
  wget https://github.com/ThePhD/sol2/releases/download/v3.2.2/sol.hpp
  wget https://github.com/ThePhD/sol2/releases/download/v3.2.2/forward.hpp
  wget https://github.com/ThePhD/sol2/releases/download/v3.2.2/config.hpp
  cd ../
fi

export GMSHSDK=${PWD}/gmsh-4.8.3-Linux64-sdk/
export EIGENSDK=${PWD}/eigen-3.3.9/
export CGAL_DIR=${PWD}/CGAL-5.2.1/
export SOLSDK=${PWD}

export PATH=${GMSHSDK}/bin:${GMSHSDK}/lib:"${PATH}"
export INCLUDE=${GMSHSDK}/include:"${INCLUDE}"
export INCLUDE=${SOLSDK}:"${INCLUDE}"
export INCLUDE=${EIGENSDK}:"${INCLUDE}"
export LIB=${GMSHSDK}/lib:"${LIB}"

cd ../

rm -rf build
mkdir build
cd build

cmake ../ -DCMAKE_BUILD_TYPE=Release -DUSE_MARCH=1  -G "Unix Makefiles"

make

cd ../