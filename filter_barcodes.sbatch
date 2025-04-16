#!/bin/bash

#### SLURM Parameters: ####
#SBATCH --nodes=1
#SBATCH --ntasks=16             # Adjust ntasks based on the Python script's ability to parallelize, or just keep 1 if it's single-threaded.
#SBATCH --job-name=extract_bc  # Changed job name slightly
#SBATCH --time=00-4:00:00
#SBATCH --mem=128GB
#SBATCH --mail-user=''          # <-- IMPORTANT: Add your email address here
#SBATCH --mail-type=ALL
#SBATCH --output='extract-barcodes-%j.out'
#SBATCH --error='extract-barcodes-%j.error'
#SBATCH --partition=prod        # Make sure 'prod' is the correct partition name

#### Script Usage ####
# This script extracts barcodes from FASTQ files using a specified pattern.
# Usage: sbatch <this_script.sh> <input_fastq> <output_prefix> [pattern_type]
#   <input_fastq>: Path to the input FASTQ file.
#   <output_prefix>: Prefix for the output files (e.g., 'sample1_barcodes'). The script will likely create files like <output_prefix>.fastq, etc.
#   [pattern_type]: Optional. Type of barcode pattern. Options: GFP, BFP, mCHERRY. Defaults to GFP.

#### Load necessary modules (if required) ####
# module purge # Uncomment if you need a clean environment
# module load python/3.x # Example: Load a specific Python module if needed

#### RUN Parameters: ####

# --- Input Arguments ---
INFILE=$1
OUTPUT=$2
# Set PATTERN_TYPE from the 3rd argument, default to 'GFP' if not provided
PATTERN_TYPE=${3:-GFP}

# --- Script Location ---
# Ensure this path is correct and the script is executable
SCRIPT="/home/YangZongmian/SPLINTR/SPLINTR_preprocessing/extractBarcodeReads.py"

# --- Define Barcode Patterns ---
PATTERN_GFP='([ATCG][ATCG][GC][AT][GC][ATGC][ATCG][AT][GC][AT]){6}'     # NNSWS NNWSW - SPLINTR GFP pattern
PATTERN_BFP='([ATCG][ATCG][AT][GC][AT][ATGC][ATCG][GC][AT][GC]){6}'     # NNWSW NNSWS - SPLINTR BFP pattern
PATTERN_MCHERRY='([ATCG][AT][ATGC][GC][ATGC][AT][ATCG][GC][ATGC][AT][ATCG][GC]){5}' # NWNSNW NSNWNS - SPLINTR mCHERRY pattern

# --- Constant Sequences ---
UPCONSTANT='TACGATTGACTA'
DOWNCONSTANT='TGCTAATGCGTACTG'

# --- Other Parameters ---
QUALITY_THRESHOLD=20
LENGTH_THRESHOLD=60 # This is the barcode length *expected* by the regex

# --- Argument Validation ---
if [ -z "$INFILE" ]; then
  echo "Error: Input FASTQ file not specified."
  echo "Usage: sbatch $0 <input_fastq> <output_prefix> [pattern_type]"
  exit 1
fi

if [ -z "$OUTPUT" ]; then
  echo "Error: Output prefix not specified."
  echo "Usage: sbatch $0 <input_fastq> <output_prefix> [pattern_type]"
  exit 1
fi

if [ ! -f "$SCRIPT" ]; then
    echo "Error: Python script not found at $SCRIPT"
    exit 1
fi

if [ ! -f "$INFILE" ]; then
    echo "Error: Input file not found at $INFILE"
    exit 1
fi


# --- Select the Pattern based on PATTERN_TYPE ---
PATTERN="" # Initialize PATTERN variable
case $PATTERN_TYPE in
  GFP)
    PATTERN=$PATTERN_GFP
    echo "Selected pattern type: GFP"
    ;;
  BFP)
    PATTERN=$PATTERN_BFP
    echo "Selected pattern type: BFP"
    ;;
  mCHERRY)
    PATTERN=$PATTERN_MCHERRY
    # Note: mCHERRY pattern length is 12*5 = 60, same as GFP/BFP (10*6 = 60)
    echo "Selected pattern type: mCHERRY"
    ;;
  *)
    echo "Error: Invalid pattern type '$PATTERN_TYPE'. Choose from GFP, BFP, mCHERRY."
    exit 1
    ;;
esac

echo "Using regex pattern: $PATTERN"

# --- Print Job Information ---
echo "Starting job $SLURM_JOB_ID on $SLURM_JOB_NODELIST"
echo "Input file: $INFILE"
echo "Output prefix: $OUTPUT"
echo "Pattern Type: $PATTERN_TYPE"
echo "Upstream Constant: $UPCONSTANT"
echo "Downstream Constant: $DOWNCONSTANT"
echo "Quality Threshold: $QUALITY_THRESHOLD"
echo "Expected Barcode Length: $LENGTH_THRESHOLD" # Check if this matches the -l param meaning in the python script
echo "Python Script: $SCRIPT"

# --- Execute Python Script ---
echo "Running Python script..."

# Ensure Python executable is found (use 'python3' if needed)
# Add full path to python if necessary, e.g., /usr/bin/python3
python "$SCRIPT" \
  --input "$INFILE" \
  --output "$OUTPUT" \
  --regex "$PATTERN" \
  --upconstant "$UPCONSTANT" \
  --downconstant "$DOWNCONSTANT" \
  -q "$QUALITY_THRESHOLD" \
  -l "$LENGTH_THRESHOLD"  # Assuming -l corresponds to the expected barcode length

EXIT_CODE=$? # Capture exit code of the python script

if [ $EXIT_CODE -eq 0 ]; then
  echo "Python script finished successfully."
else
  echo "Python script failed with exit code $EXIT_CODE."
  # Optional: Add commands here if you need to do something specific upon failure
fi

echo "Job finished."

exit $EXIT_CODE