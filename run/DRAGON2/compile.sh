#!/bin/bash
#Do not forget to install libgfortran3  on your system!

module load CMake/3.12.1-GCCcore-7.3.0
module load CGAL/4.11.1-foss-2018b-Python-2.7.15

cd ../../

if [ ! -d "gmsh-4.8.3-Linux64-sdk" ]; then
  wget http://gmsh.info/bin/Linux/gmsh-4.8.3-Linux64-sdk.tgz
  tar -xf gmsh-4.8.3-Linux64-sdk.tgz 
  rm -rf gmsh-4.8.3-Linux64-sdk.tgz 
fi

if [ ! -d "eigen-eigen-323c052e1731" ]; then
  wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz
  tar -xf 3.3.7.tar.gz
  rm -rf 3.3.7.tar.gz
fi

if [ ! -d "nlohmann" ]; then
  wget https://github.com/nlohmann/json/releases/download/v3.6.1/include.zip
  unzip include.zip
  rm -rf include.zip
  mkdir nlohmann/
  mv include/nlohmann nlohmann/nlohmann
  rm -rf include
fi

export GMSHSDK=${PWD}/gmsh-4.8.3-Linux64-sdk/
export EIGENSDK=${PWD}/eigen-eigen-323c052e1731/
export JSONSDK=${PWD}/nlohmann/:

export PATH=${GMSHSDK}/bin:${GMSHSDK}/lib:"${PATH}"
export INCLUDE=${GMSHSDK}/include:"${INCLUDE}"
export INCLUDE=${EIGENSDK}:"${INCLUDE}"
export INCLUDE=${JSONSDK}:"${INCLUDE}"
export LIB=${GMSHSDK}/lib:"${LIB}"
export PYTHONPATH=${GMSHSDK}/lib:"${PYTHONPATH}"
export DYLD_LIBRARY_PATH=${GMSHSDK}/lib:"${DYLD_LIBRARY_PATH}"

rm -rf build
mkdir build
cd build

cmake ../ -DCMAKE_BUILD_TYPE=Release -DUSE_MARCH=1  -G "Unix Makefiles"

make

cd ../