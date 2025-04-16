#!/bin/bash

#### SLURM Parameters: ####
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --job-name=flash
#SBATCH --time=00-10:00:00
#SBATCH --mem=256GB
#SBATCH --mail-user=''
#SBATCH --mail-type=ALL
#SBATCH --output='flash-%j.out'
#SBATCH --error='flash-%j.error'
#SBATCH --partition=prod

#### Load default modules ####
module purge
module load java

#### RUN Parameters: ####
READ1=$1
READ2=$2
OUTDIR="../flash"

# get output folder name and create if req
name=$(basename $READ1)
name=${name%%.*}
name=$(echo $name | cut -d "_" -f1)
echo $name

OUT=${OUTDIR}/${name}

if [ ! -d $OUT ]; then
    mkdir $OUT
fi

# place full path of flash software here 
FLASH="place full path to flash executable here"

# print version of flash
$FLASH -v

# run flash with defaults on read1 and read2 
cmd="$FLASH -d $OUT -o $name -z -t 16 $READ1 $READ2"
echo $cmd
$cmd


