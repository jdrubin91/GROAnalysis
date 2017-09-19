###Name the job
#PBS -N GSEA
### Specify the number of nodes/cores
#PBS -l nodes=1:ppn=1

### Allocate the amount of memory needed
#PBS -l mem=8gb

### Set your expected walltime
#PBS -l walltime=05:00:00

#PBS -q short8gb

### Setting to mail when the job is complete
#PBS -e /Users/joru1876/qsub_errors/
#PBS -o /Users/joru1876/qsub_stdo/  

### Set your email address
#PBS -m ae
#PBS -M joru1876@colorado.edu



### Choose your shell 
#PBS -S /bin/bash
### Pass enviroment variables to the job

module load bedtools2_2.22.0

### now call your program

src=/Users/joru1876/scratch_backup/GROAnalysis/src/GSEA.py

python $src

