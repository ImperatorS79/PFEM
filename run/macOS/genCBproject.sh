#!/bin/sh
#Do not forget to install libgfortran3  on your system!
#brew update && brew install cmake wget gnu-tar llvm libomp unzip swig cgal lua python@3.7 nlohmann-json eigen
#This project needs a C++17 compliant compiler

cd ../../
if [ ! -d "dependencies" ]; then
    mkdir dependencies
fi

cd dependencies/

if [ ! -d "sol" ]; then
  mkdir sol
  cd sol
  wget https://github.com/ThePhD/sol2/releases/download/v3.2.2/sol.hpp
  wget https://github.com/ThePhD/sol2/releases/download/v3.2.2/forward.hpp
  wget https://github.com/ThePhD/sol2/releases/download/v3.2.2/config.hpp
  cd ../
fi

if [ ! -d "gmsh-4.7.1-MacOSX-sdk.tgz" ]; then
  wget http://gmsh.info/bin/MacOSX/gmsh-4.7.1-MacOSX-sdk.tgz
  gtar -xf gmsh-4.7.1-MacOSX-sdk.tgz 
  rm -rf gmsh-4.7.1-MacOSX-sdk.tgz 
fi

export GMSHSDK=${PWD}/gmsh-4.7.1-MacOSX-sdk/
export SOLSDK=${PWD}/
export EIGENSDK=/usr/local/include/eigen3/


export PATH=${GMSHSDK}/bin:${GMSHSDK}/lib:"${PATH}"
export INCLUDE=${GMSHSDK}/include:"${INCLUDE}"
export INCLUDE=${SOLSDK}:"${INCLUDE}"
export INCLUDE=${EIGENSDK}:"${INCLUDE}"
export LIB=${GMSHSDK}/lib:"${LIB}"
export PYTHONPATH=${GMSHSDK}/lib:"${PYTHONPATH}"
export DYLD_LIBRARY_PATH=${GMSHSDK}/lib:"${DYLD_LIBRARY_PATH}"

cd ../

rm -rf build
mkdir build
cd build

LDFLAGS="-L/usr/local/opt/llvm/lib -L/usr/local/opt/lua/lib/" CPPFLAGS="-I/usr/local/opt/llvm/include -I/usr/local/opt/lua/include/" cmake ../ -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=/usr/local/opt/llvm/bin/clang -DCMAKE_CXX_COMPILER=/usr/local/opt/llvm/bin/clang++  -G "CodeBlocks - Unix Makefiles"

cp ../run/macos/run.sh bin/

cd ../
