#!/bin/bash

#### SLURM Parameters: ####
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --job-name=fastqc
#SBATCH --time=00-04:00:00
#SBATCH --mem=64GB
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH --output='fastqc-%j.out'
#SBATCH --error='fastqc-%j.error'

#### Load default modules ####
module purge
module load java
module load fastqc

#### RUN Parameters: ####
INPUT=$1
fastqc $INPUT -t 8 -o "../qc"
