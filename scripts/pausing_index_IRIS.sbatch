#!/bin/bash
#SBATCH -p long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output /scratch/Users/joru1876/GROAnalysis/e_and_o/%x.out
#SBATCH --error /scratch/Users/joru1876/GROAnalysis/e_and_o/%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=joru1876@colorado.edu
module load bedtools/2.25.0
module load python/2.7.14
module load python/2.7.14/matplotlib/1.5.1
module load python/2.7.14/numpy/1.11.0
module load python/2.7.14/scipy/0.17.1

src=/scratch/Users/joru1876/GROAnalysis/src/pausing_index_IRIS.py

python $src
