#!/bin/bash
#SBATCH -n 1
#SBATCH -t 00:00:15
#SBATCH --mem-per-cpu=2500
#SBATCH -p debug

echo "Starting sphere decoder simulation..."
echo "Loading module Armadillo/7.800.2-iomkl-triton-2017a-Python-2.7.13"
module restore arma-env

# g++ sphdec.cpp -o sphdec -O2 -larmadillo -llapack -lblas -std=c++14
cd ..
make
srun pwalk
