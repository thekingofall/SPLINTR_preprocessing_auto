#!/bin/bash

#### SLURM Parameters: ####
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --job-name=starcode
#SBATCH --time=00-04:00:00
#SBATCH --mem=256GB
#SBATCH --mail-user=''
#SBATCH --mail-type=ALL
#SBATCH --output='starcode-%j.out'
#SBATCH --error='starcode-%j.error'
#SBATCH --partition=prod

#### Load default modules ####
# module purge
# module load java
# module load fastqc

#### RUN Parameters: ####
INPUT=$1
OUTDIR="../starcode"

name=$(basename $INPUT)
name=${name%%.*}


STARCODE="full path to starcode executable goes here"
$STARCODE -v

# distance auto
$STARCODE -t 8 $INPUT > ${OUTDIR}/${name}_starcode_auto.txt
