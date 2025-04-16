#!/bin/bash


# if using indexes place full path to that text file here
INDEXES=""

# place full path to desired output directory here
OUTDIR="./" 

# setup output directory if it does not already exist
if [ ! -d $OUTDIR ]; then
    mkdir $OUTDIR
fi

# lauch extractBarcodeReads script onto SLURM HPC

# change following find command to path containing fastq files 
for file in $(find ../flash -name '*extendedFrags.fastq.gz' -type 'f')
    do
        outname=$(basename $file)
        outname=${outname%%.*}_filtered.fastq
        dir=$(dirname $file)
        #echo $file
        #echo $outname
        #echo $dir
        cmd="sbatch filter_barcodes.sbatch $file ${OUTDIR}/${outname}"
        echo $cmd

        # uncomment the below line once you validate that the command works as expected.
        $cmd
    done

