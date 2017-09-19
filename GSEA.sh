###Name the job
#PBS -N name your job
### Specify the number of nodes/cores
#PBS -l nodes=1:ppn=1

### Allocate the amount of memory needed
#PBS -l mem=10gb

### Set your expected walltime
#PBS -l walltime=100:00:00

### Setting to mail when the job is complete
#PBS -e /scratch/Users/joru1876/qsub_errors/
#PBS -o /scratch/Users/joru1876/qsub_stdo/  

### Set your email address
#PBS -m ae
#PBS -M joru1876@colorado.edu



### Choose your shell 
#PBS -S /bin/bash
### Pass enviroment variables to the job

module load bedtools2_2.22.0

### now call your program

src=/Users/joru1876/scratch_backup/GROAnalysis/GSEA.py

python $src

