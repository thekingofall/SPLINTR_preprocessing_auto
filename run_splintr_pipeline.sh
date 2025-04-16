#!/bin/bash

# --- Default Settings ---
DEFAULT_OUTPUT_DIR="./Processed_SPLINTR_Output_Parallel"
STARCODE_THREADS=2 # Default threads *per starcode instance* when run in parallel
CONDA_ENV="SPLINTR" # The name of the conda environment to use
DEFAULT_PARALLEL_JOBS=12 # Default number of samples to process concurrently
DEFAULT_RUN_ANALYSIS="true" # Default is to run the R analysis script
# IMPORTANT: Make sure these paths are correct for your system
SBATCH_SCRIPT="$HOME/SPLINTR/SPLINTR_preprocessing/filter_barcodes.sbatch"
# IMPORTANT: Make sure this path is correct for your R analysis script
R_ANALYSIS_SCRIPT="$HOME/SPLINTR/SPLINTR_preprocessing/run_splintr_draw.R"


# --- Functions ---
print_usage() {
  echo "Usage: $0 -i <input_directory> [-o <output_directory>] [-t <starcode_threads>] [-p <parallel_jobs>] [-r <true|false>] [-h]"
  echo ""
  echo "Processes SPLINTR R1 FASTQ files in parallel, runs starcode, and optionally runs a final R analysis script."
  echo ""
  echo "Options:"
  echo "  -i <input_directory>    : Directory containing sample subdirectories (e.g., 1/, 2/, 10/)."
  echo "                            Each subdirectory must contain a corresponding _1.fq.gz file."
  echo "                            (Mandatory)"
  echo "  -o <output_directory>   : Main directory where processed files and logs will be saved."
  echo "                            Step results will be in sub-folders (Step1_Barcodes, Step2_Starcode)."
  echo "                            (Optional, default: ${DEFAULT_OUTPUT_DIR})"
  echo "  -t <starcode_threads> : Number of threads for each starcode instance."
  echo "                            (Optional, default: ${STARCODE_THREADS})"
  echo "  -p <parallel_jobs>    : Number of samples to process in parallel via ParaFly."
  echo "                            (Optional, default: ${DEFAULT_PARALLEL_JOBS})"
  echo "  -r <true|false>       : Run the final R analysis script (run_splintr_draw.R) on the Step2 output."
  echo "                            (Optional, default: ${DEFAULT_RUN_ANALYSIS})"
  echo "  -h                      : Display this help message."
  echo ""
  echo "Example:"
  echo "  # Process and run analysis (default)"
  echo "  $0 -i /path/to/01.RawData -o /path/to/results -p 8 -t 2"
  echo ""
  echo "  # Process but DO NOT run analysis"
  echo "  $0 -i /path/to/01.RawData -p 10 -r false"
}

# --- Argument Parsing ---
INPUT_DIR=""
OUTPUT_DIR="$DEFAULT_OUTPUT_DIR" # Set default
PARALLEL_JOBS="$DEFAULT_PARALLEL_JOBS"
RUN_ANALYSIS="$DEFAULT_RUN_ANALYSIS" # Set default

while getopts ":i:o:t:p:r:h" opt; do
  case ${opt} in
    i )
      INPUT_DIR=$OPTARG
      ;;
    o )
      OUTPUT_DIR=$OPTARG
      ;;
    t )
      STARCODE_THREADS=$OPTARG
      # Validate threads input immediately
      if ! [[ "$STARCODE_THREADS" =~ ^[0-9]+$ ]] || [ "$STARCODE_THREADS" -lt 1 ]; then
          echo "Error: Starcode threads (-t) must be a positive integer. Found: '$STARCODE_THREADS'"
          print_usage
          exit 1
      fi
      ;;
    p )
      PARALLEL_JOBS=$OPTARG
       # Validate parallel jobs input immediately
      if ! [[ "$PARALLEL_JOBS" =~ ^[0-9]+$ ]] || [ "$PARALLEL_JOBS" -lt 1 ]; then
          echo "Error: Parallel jobs (-p) must be a positive integer. Found: '$PARALLEL_JOBS'"
          print_usage
          exit 1
      fi
      ;;
    r )
      RUN_ANALYSIS=$(echo "$OPTARG" | tr '[:upper:]' '[:lower:]') # Convert to lowercase
      if [[ "$RUN_ANALYSIS" != "true" && "$RUN_ANALYSIS" != "false" ]]; then
        echo "Error: Invalid value for -r option. Must be 'true' or 'false'. Found: '$OPTARG'"
        print_usage
        exit 1
      fi
      ;;
    h )
      print_usage
      exit 0
      ;;
    \? )
      echo "Invalid option: $OPTARG" 1>&2
      print_usage
      exit 1
      ;;
    : )
      echo "Invalid option: $OPTARG requires an argument" 1>&2
      print_usage
      exit 1
      ;;
  esac
done
shift $((OPTIND -1))

# --- Validate Input ---
if [ -z "$INPUT_DIR" ]; then
  echo "Error: Input directory (-i) is mandatory."
  print_usage
  exit 1
fi

# Convert to absolute paths early
INPUT_DIR_ABS=$(cd "$INPUT_DIR" && pwd)
if [ $? -ne 0 ] || [ ! -d "$INPUT_DIR_ABS" ]; then
  echo "Error: Input directory '$INPUT_DIR' not found or is not a directory."
  exit 1
fi
INPUT_DIR="$INPUT_DIR_ABS" # Use absolute path

# Final check on numeric defaults if not overridden
if ! [[ "$STARCODE_THREADS" =~ ^[0-9]+$ ]] || [ "$STARCODE_THREADS" -lt 1 ]; then
    echo "Error: Default starcode threads value is invalid. Check script default settings."
    exit 1
fi
if ! [[ "$PARALLEL_JOBS" =~ ^[0-9]+$ ]] || [ "$PARALLEL_JOBS" -lt 1 ]; then
    echo "Error: Default parallel jobs value is invalid. Check script default settings."
    exit 1
fi

# --- Prepare Output Directories ---
# Make OUTPUT_DIR absolute *before* creating it and subdirs
OUTPUT_DIR_ABS=$(mkdir -p "$OUTPUT_DIR" && cd "$OUTPUT_DIR" && pwd)
if [ $? -ne 0 ] || [ -z "$OUTPUT_DIR_ABS" ]; then
  echo "Error: Failed to create or access main output directory '$OUTPUT_DIR'."
  exit 1
fi
OUTPUT_DIR="$OUTPUT_DIR_ABS" # Use absolute path henceforth

# Define and create step-specific sub-directories
STEP1_OUT_DIR="${OUTPUT_DIR}/Step1_Barcodes"
STEP2_OUT_DIR="${OUTPUT_DIR}/Step2_Starcode"
# Define analysis output directory (used by the R script)
ANALYSIS_OUT_DIR="${OUTPUT_DIR}/Step3_SPLINTR_Analysis_Output"

echo "INFO: Creating sub-directories for steps:"
echo "      Step 1: $STEP1_OUT_DIR"
echo "      Step 2: $STEP2_OUT_DIR"
echo "      Step 3 (Analysis): $ANALYSIS_OUT_DIR" # Also list analysis dir
mkdir -p "$STEP1_OUT_DIR" "$STEP2_OUT_DIR" "$ANALYSIS_OUT_DIR" # Create all needed dirs
if [ $? -ne 0 ]; then
    echo "Error: Failed to create step sub-directories in '$OUTPUT_DIR'."
    exit 1
fi

# --- Setup Logging ---
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
# Place log file in the main output directory
LOG_FILE="${OUTPUT_DIR}/splintr_pipeline_${TIMESTAMP}.log"
echo "INFO: Script execution log will be saved to: $LOG_FILE"
# Redirect all subsequent stdout and stderr to tee, appending to the log file
exec > >(tee -a "$LOG_FILE") 2> >(tee -a "$LOG_FILE" >&2)

echo "INFO: Logging started at $(date)"
echo "=================================================="

# --- Verify Commands and Scripts ---
echo "INFO: Verifying required commands and scripts..."
COMMAND_NOT_FOUND=0
if ! command -v conda &> /dev/null; then echo "Error: 'conda' not found."; COMMAND_NOT_FOUND=1; fi
if ! command -v ParaFly &> /dev/null; then echo "Error: 'ParaFly' not found."; COMMAND_NOT_FOUND=1; fi
if ! command -v Rscript &> /dev/null; then echo "Error: 'Rscript' not found."; COMMAND_NOT_FOUND=1; fi

# Ensure SBATCH_SCRIPT and R_ANALYSIS_SCRIPT are absolute paths
if [[ "$SBATCH_SCRIPT" != /* ]]; then SBATCH_SCRIPT_ABS=$(cd "$(dirname "$SBATCH_SCRIPT")" && pwd)/$(basename "$SBATCH_SCRIPT"); else SBATCH_SCRIPT_ABS="$SBATCH_SCRIPT"; fi
if [ ! -f "$SBATCH_SCRIPT_ABS" ]; then echo "Error: SBATCH script not found: $SBATCH_SCRIPT"; COMMAND_NOT_FOUND=1; else SBATCH_SCRIPT="$SBATCH_SCRIPT_ABS"; fi

if [[ "$R_ANALYSIS_SCRIPT" != /* ]]; then R_ANALYSIS_SCRIPT_ABS=$(cd "$(dirname "$R_ANALYSIS_SCRIPT")" && pwd)/$(basename "$R_ANALYSIS_SCRIPT"); else R_ANALYSIS_SCRIPT_ABS="$R_ANALYSIS_SCRIPT"; fi
if [ ! -f "$R_ANALYSIS_SCRIPT_ABS" ]; then echo "Error: R analysis script not found: $R_ANALYSIS_SCRIPT"; COMMAND_NOT_FOUND=1; else R_ANALYSIS_SCRIPT="$R_ANALYSIS_SCRIPT_ABS"; fi

if [ $COMMAND_NOT_FOUND -eq 1 ]; then
    echo "Error: One or more required commands or scripts were not found. Please check paths and installations."
    exit 1
fi
echo "INFO: All required commands and scripts found."

# --- Prepare Command List ---
COMMAND_FILE=$(mktemp "${OUTPUT_DIR}/splintr_parafly_cmds.XXXXXX")
if [ $? -ne 0 ] || [ -z "$COMMAND_FILE" ] || [ ! -f "$COMMAND_FILE" ]; then
    echo "Error: Failed to create temporary command file in $OUTPUT_DIR"
    exit 1
fi
echo "INFO: Generating commands list in: $COMMAND_FILE"

num_commands=0

echo "=================================================="
echo "Configuration:"
echo "Input Directory: $INPUT_DIR"
echo "Output Directory (Main): $OUTPUT_DIR"
echo "  -> Barcodes Output: $STEP1_OUT_DIR"
echo "  -> Starcode Output: $STEP2_OUT_DIR"
echo "  -> Analysis Output: $ANALYSIS_OUT_DIR"
echo "Starcode Threads per job: $STARCODE_THREADS"
echo "Parallel Jobs (ParaFly -CPU): $PARALLEL_JOBS"
echo "SBATCH Script: $SBATCH_SCRIPT"
echo "Conda Environment: $CONDA_ENV"
echo "Run Final Analysis: $RUN_ANALYSIS"
if [ "$RUN_ANALYSIS" = "true" ]; then
    echo "R Analysis Script: $R_ANALYSIS_SCRIPT"
fi
echo "Command File: $COMMAND_FILE"
echo "Log File: $LOG_FILE"
echo "=================================================="
echo "Generating commands for Steps 1 & 2..."

shopt -s nullglob # Prevent loop from running if no matches found

for sample_dir_path in "$INPUT_DIR"/*/; do
    sample_name=$(basename "$sample_dir_path")
    # Check if the directory name looks like a sample identifier
    if [[ "$sample_name" =~ ^[a-zA-Z0-9_-]+$ ]]; then
        input_fastq_abs="${INPUT_DIR}/${sample_name}/${sample_name}_1.fq.gz"
        if [ -f "$input_fastq_abs" ]; then
            step1_output_abs="${STEP1_OUT_DIR}/Sample${sample_name}_barcodes.txt"
            step2_output_abs="${STEP2_OUT_DIR}/Sample${sample_name}_starcode.txt"
            cmd_str="echo '--- Starting Sample: ${sample_name} ---' && "
            cmd_str+="( conda run -n \"${CONDA_ENV}\" --no-capture-output --live-stream bash \"${SBATCH_SCRIPT}\" \"${input_fastq_abs}\" \"${step1_output_abs}\" ) && "
            cmd_str+="( conda run -n \"${CONDA_ENV}\" --no-capture-output --live-stream starcode -t ${STARCODE_THREADS} \"${step1_output_abs}\" > \"${step2_output_abs}\" ) && "
            cmd_str+="echo '--- Finished Sample: ${sample_name} ---' || "
            cmd_str+="echo '*** Error processing Sample: ${sample_name} ***'"
            echo "$cmd_str" >> "$COMMAND_FILE"
            num_commands=$((num_commands + 1))
        else
            echo "Warning: Input file $input_fastq_abs not found. Skipping sample $sample_name."
        fi
    else
         echo "Skipping directory (doesn't look like a sample): $sample_dir_path"
    fi
done
shopt -u nullglob # Turn off nullglob

if [ $num_commands -eq 0 ]; then
    echo "Error: No valid sample directories containing '_1.fq.gz' files found in $INPUT_DIR"
    rm -f "$COMMAND_FILE" # Clean up empty temp file
    exit 1
fi

echo "INFO: Generated $num_commands commands for ParaFly."
echo "=================================================="
echo "Starting parallel processing (Steps 1 & 2) with ParaFly at $(date)..."
echo "Output and errors from ParaFly will follow (also logged to $LOG_FILE):"
echo "=================================================="

# --- Execute Steps 1 & 2 with ParaFly ---
FAILED_CMDS_LOG="${OUTPUT_DIR}/parafly_failed_cmds_${TIMESTAMP}.log"
ParaFly -c "$COMMAND_FILE" -CPU "$PARALLEL_JOBS" -vv -failed_cmds "$FAILED_CMDS_LOG"
parafly_exit_code=$?

echo "=================================================="
echo "ParaFly execution (Steps 1 & 2) finished at $(date)."

# --- Check ParaFly Results ---
PARAFFLY_FAILED=0
if [ $parafly_exit_code -ne 0 ]; then
    echo "Error: ParaFly reported an error (Exit code: $parafly_exit_code)."
    PARAFFLY_FAILED=1
    if [ -s "$FAILED_CMDS_LOG" ]; then
        echo "       Check the failed commands log: $FAILED_CMDS_LOG"
    else
        echo "       No failed commands log was generated or it is empty."
    fi
else
    echo "ParaFly completed successfully for Steps 1 & 2."
    if [ -s "$FAILED_CMDS_LOG" ]; then
        echo "Warning: ParaFly exited successfully, but the failed commands log is not empty:"
        echo "         $FAILED_CMDS_LOG"
        echo "         This might indicate issues within successfully completed ParaFly jobs (check specific command outputs)."
        # Decide if this warning should prevent analysis - currently it doesn't
    else
         # rm -f "$FAILED_CMDS_LOG" # Optionally remove the empty failed cmds log on success
         : # Bash no-op if you don't want to remove it
    fi
fi

# --- Cleanup ParaFly Command File ---
echo "INFO: Removing temporary command file: $COMMAND_FILE"
rm -f "$COMMAND_FILE"

# --- Optionally Run R Analysis Script (Step 3) ---
FINAL_EXIT_CODE=$parafly_exit_code # Start with ParaFly's exit code

if [ "$RUN_ANALYSIS" = "true" ]; then
    echo "=================================================="
    echo "Proceeding to run R analysis script (Step 3)..."

    # Check if ParaFly failed - perhaps don't run R if it did? (Optional decision)
    if [ $PARAFFLY_FAILED -eq 1 ]; then
        echo "Warning: ParaFly reported errors in Steps 1/2. Analysis in Step 3 might be based on incomplete data."
        # You could add 'exit $parafly_exit_code' here if you want to stop
    fi

    # Check if the input directory for R script (Step2 output) exists and is not empty
    if [ ! -d "$STEP2_OUT_DIR" ] || [ -z "$(ls -A "$STEP2_OUT_DIR"/*.txt 2>/dev/null)" ]; then
         echo "Error: Step 2 output directory '$STEP2_OUT_DIR' is missing or empty. Cannot run R analysis."
         # Update exit code if not already failed
         if [ $FINAL_EXIT_CODE -eq 0 ]; then FINAL_EXIT_CODE=1; fi
    else
        echo "INFO: Running R analysis script:"
        echo "      Input: $STEP2_OUT_DIR"
        echo "      Output: $ANALYSIS_OUT_DIR"
        echo "      Log: Appended to $LOG_FILE"
        echo "--------------------------------------------------"

        # Use absolute paths for R script arguments
        # Note: The R script itself handles the specific log file name within its output dir
        # Ensure the R script is executable or use 'Rscript' command
        Rscript "$R_ANALYSIS_SCRIPT" --input_dir "$STEP2_OUT_DIR" --output_dir "$ANALYSIS_OUT_DIR" # Add other R script args if needed
        r_exit_code=$?
        echo "--------------------------------------------------"
        echo "R analysis script finished at $(date) with exit code: $r_exit_code"

        if [ $r_exit_code -ne 0 ]; then
            echo "Error: R analysis script failed (Exit code: $r_exit_code)."
            # Update overall exit code if the script succeeded up to this point
            if [ $FINAL_EXIT_CODE -eq 0 ]; then FINAL_EXIT_CODE=$r_exit_code; fi
        else
            echo "R analysis script completed successfully."
        fi
    fi
else
    echo "=================================================="
    echo "INFO: Skipping final R analysis step as requested (-r false)."
fi


# --- Final Summary ---
echo "=================================================="
echo "Processing Complete at $(date)."
echo "Results saved in main directory: $OUTPUT_DIR"
echo " -> Step 1 (Barcodes) results in: $STEP1_OUT_DIR"
echo " -> Step 2 (Starcode) results in: $STEP2_OUT_DIR"
if [ "$RUN_ANALYSIS" = "true" ]; then
    echo " -> Step 3 (Analysis) results in: $ANALYSIS_OUT_DIR"
fi
echo "Main execution log: $LOG_FILE"
if [ -f "$FAILED_CMDS_LOG" ]; then # Check if failed log exists before mentioning it
    echo "Failed ParaFly commands log (if any): $FAILED_CMDS_LOG"
fi
echo "=================================================="

# Exit with the final exit code (either from ParaFly or R script if it failed later)
exit $FINAL_EXIT_CODE