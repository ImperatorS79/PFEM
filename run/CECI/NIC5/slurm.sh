#!/bin/bash
# Submission script for Nic5
#SBATCH --job-name=DropFallInFluid
#SBATCH --time=48:00:00 # hh:mm:ss
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=10 # megabytes
#SBATCH --partition=batch
#
#SBATCH --mail-user=simon.fevrier@uliege.be
#SBATCH --mail-type=ALL
#
#SBATCH --comment=PFEM3D
#
#SBATCH --output=DropFallInFluid.txt

export OMP_NUM_THREADS=64
export MKL_NUM_THREADS=64

export PATH=../../dependencies/gmsh-4.8.3-Linux64-sdk/bin:"${PATH}"

./pfem examples/3D/dropFallInFluid/testComp.lua