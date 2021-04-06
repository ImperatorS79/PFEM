#!/bin/sh
#Launch this in MinGW64 shell from MSYS2

buildType="Release"
if [ ! -z "$2" ]; then
	buildType="$2"
fi

cd ../../
if [ ! -d "dependencies" ]; then
    mkdir dependencies
fi

cd dependencies/

if [ ! -d "gmsh-4.7.0-Windows64-sdk" ]; then
  wget http://gmsh.info/bin/Windows/gmsh-4.7.0-Windows64-sdk.zip
  unzip gmsh-4.7.0-Windows64-sdk.zip 
  rm -rf gmsh-4.7.0-Windows64-sdk.zip
fi

export CGAL_DIR=/c/tools/msys64/mingw64/lib/cmake/CGAL/
export GMSHSDK=${PWD}/gmsh-4.7.0-Windows64-sdk #put gmsh sdk here
export EIGENSDK=/c/tools/msys64/mingw64/include/eigen3
export SOL3SDK=/c/tools/mingw64/include/

export PATH=${GMSHSDK}/bin:${GMSHSDK}/lin:${PATH}
export INCLUDE=${GMSHSDK}/include:${EIGENSDK}:${SOL3SDK}
export LIB=${GMSHSDK}/lib

cd ../

rm -rf build
mkdir build
cd build

cmake -G "CodeBlocks - MinGW Makefiles" -DCGAL_HEADER_ONLY=ON -DUSE_MARCH=0 -DCMAKE_SH=SH-NOTFOUND  ..

cd ../
