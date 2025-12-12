#!/bin/bash

# Script to test DANTE_TIR with different parameter combinations
# Usage: ./run_parameter_tests.sh <test-data-dir> <output-dir>

set -u  # Exit on undefined variables

# Check arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 <test-data-dir> <output-dir>"
    echo "Example: $0 test-data output"
    echo ""
    echo "Arguments:"
    echo "  test-data-dir  Directory containing Ath, Med, Ghi, Sce subdirectories"
    echo "  output-dir     Directory where results will be stored"
    exit 1
fi

TEST_DATA_DIR="$1"
OUTPUT_DIR="$2"

# Check if test data directory exists
if [ ! -d "$TEST_DATA_DIR" ]; then
    echo "Error: Test data directory '$TEST_DATA_DIR' does not exist"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Parameters to test
GENOMES=("Ath" "Med" "Ghi" "Sce")
N_BEAST_ITER=(1 5 10)
MAX_CLASS_SIZE=("not_used" 2000 4000)
THREADS=20

# Get the directory where dante_tir.py is located (same as this script)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DANTE_TIR="$SCRIPT_DIR/dante_tir.py"

# Check if dante_tir.py exists
if [ ! -f "$DANTE_TIR" ]; then
    echo "Error: dante_tir.py not found at $DANTE_TIR"
    exit 1
fi

# Log files
MAIN_LOG="$OUTPUT_DIR/test_run.log"
TIMING_CSV="$OUTPUT_DIR/timing_results.csv"

# Initialize main log
echo "========================================" | tee "$MAIN_LOG"
echo "DANTE_TIR Parameter Testing" | tee -a "$MAIN_LOG"
echo "========================================" | tee -a "$MAIN_LOG"
echo "Started at: $(date)" | tee -a "$MAIN_LOG"
echo "Test data directory: $TEST_DATA_DIR" | tee -a "$MAIN_LOG"
echo "Output directory: $OUTPUT_DIR" | tee -a "$MAIN_LOG"
echo "Threads per run: $THREADS" | tee -a "$MAIN_LOG"
echo "Parameters:" | tee -a "$MAIN_LOG"
echo "  n_beast_iter: ${N_BEAST_ITER[*]}" | tee -a "$MAIN_LOG"
echo "  max_class_size: ${MAX_CLASS_SIZE[*]}" | tee -a "$MAIN_LOG"
echo "========================================" | tee -a "$MAIN_LOG"
echo "" | tee -a "$MAIN_LOG"

# Initialize timing CSV
echo "run_number,genome,n_beast_iter,max_class_size,status,elapsed_seconds,elapsed_minutes,elapsed_formatted,start_time,end_time" > "$TIMING_CSV"

# Statistics
TOTAL_RUNS=0
SUCCESSFUL_RUNS=0
FAILED_RUNS=0
SKIPPED_RUNS=0

# Loop through genomes sequentially
for genome in "${GENOMES[@]}"; do
    echo "========================================" | tee -a "$MAIN_LOG"
    echo "Processing genome: $genome" | tee -a "$MAIN_LOG"
    echo "========================================" | tee -a "$MAIN_LOG"

    GENOME_DIR="$TEST_DATA_DIR/$genome"
    GENOME_FASTA="$GENOME_DIR/genome.fasta"
    DANTE_GFF="$GENOME_DIR/DANTE_filtered.gff3"

    # Check if genome files exist
    if [ ! -f "$GENOME_FASTA" ]; then
        echo "WARNING: $GENOME_FASTA not found, skipping $genome" | tee -a "$MAIN_LOG"
        echo "" | tee -a "$MAIN_LOG"
        continue
    fi
    if [ ! -f "$DANTE_GFF" ]; then
        echo "WARNING: $DANTE_GFF not found, skipping $genome" | tee -a "$MAIN_LOG"
        echo "" | tee -a "$MAIN_LOG"
        continue
    fi

    # Loop through parameter combinations
    for n_beast in "${N_BEAST_ITER[@]}"; do
        for max_class in "${MAX_CLASS_SIZE[@]}"; do
            TOTAL_RUNS=$((TOTAL_RUNS + 1))

            # Create output directory name
            RUN_DIR="$OUTPUT_DIR/$genome/nbeast_${n_beast}_maxclass_${max_class}"
            mkdir -p "$RUN_DIR"

            # Check if output already exists (skip if completed)
            FINAL_OUTPUT="$RUN_DIR/DANTE_TIR_final.gff3"
            if [ -f "$FINAL_OUTPUT" ]; then
                SKIPPED_RUNS=$((SKIPPED_RUNS + 1))
                echo "----------------------------------------" | tee -a "$MAIN_LOG"
                echo "Run #$TOTAL_RUNS" | tee -a "$MAIN_LOG"
                echo "Genome: $genome" | tee -a "$MAIN_LOG"
                echo "Parameters: n_beast_iter=$n_beast, max_class_size=$max_class" | tee -a "$MAIN_LOG"
                echo "Output: $RUN_DIR" | tee -a "$MAIN_LOG"
                echo "STATUS: SKIPPED (output already exists)" | tee -a "$MAIN_LOG"
                echo "" | tee -a "$MAIN_LOG"
                continue
            fi

            # Log files for this run
            STDOUT_LOG="$RUN_DIR/stdout.log"
            STDERR_LOG="$RUN_DIR/stderr.log"

            # Build command
            CMD="$DANTE_TIR -g $DANTE_GFF -f $GENOME_FASTA -o $RUN_DIR -c $THREADS --n_beast_iter $n_beast"

            # Add max_class_size parameter if not "not_used"
            if [ "$max_class" != "not_used" ]; then
                CMD="$CMD --max_class_size $max_class"
            fi

            echo "----------------------------------------" | tee -a "$MAIN_LOG"
            echo "Run #$TOTAL_RUNS" | tee -a "$MAIN_LOG"
            echo "Genome: $genome" | tee -a "$MAIN_LOG"
            echo "Parameters: n_beast_iter=$n_beast, max_class_size=$max_class" | tee -a "$MAIN_LOG"
            echo "Output: $RUN_DIR" | tee -a "$MAIN_LOG"
            echo "Command: $CMD" | tee -a "$MAIN_LOG"

            # Run the command
            START_TIME=$(date +%s)
            START_TIME_FORMATTED=$(date '+%Y-%m-%d %H:%M:%S')
            echo "Started at: $START_TIME_FORMATTED" | tee -a "$MAIN_LOG"

            if $CMD > "$STDOUT_LOG" 2> "$STDERR_LOG"; then
                END_TIME=$(date +%s)
                END_TIME_FORMATTED=$(date '+%Y-%m-%d %H:%M:%S')
                ELAPSED=$((END_TIME - START_TIME))
                ELAPSED_MIN=$(echo "scale=2; $ELAPSED / 60" | bc)
                ELAPSED_FORMATTED="$(($ELAPSED / 60))m $(($ELAPSED % 60))s"
                SUCCESSFUL_RUNS=$((SUCCESSFUL_RUNS + 1))
                STATUS="SUCCESS"
                echo "STATUS: $STATUS" | tee -a "$MAIN_LOG"
                echo "Completed in: ${ELAPSED}s ($ELAPSED_FORMATTED)" | tee -a "$MAIN_LOG"

                # Record timing to CSV
                echo "$TOTAL_RUNS,$genome,$n_beast,$max_class,$STATUS,$ELAPSED,$ELAPSED_MIN,$ELAPSED_FORMATTED,$START_TIME_FORMATTED,$END_TIME_FORMATTED" >> "$TIMING_CSV"
            else
                EXIT_CODE=$?
                END_TIME=$(date +%s)
                END_TIME_FORMATTED=$(date '+%Y-%m-%d %H:%M:%S')
                ELAPSED=$((END_TIME - START_TIME))
                ELAPSED_MIN=$(echo "scale=2; $ELAPSED / 60" | bc)
                ELAPSED_FORMATTED="$(($ELAPSED / 60))m $(($ELAPSED % 60))s"
                FAILED_RUNS=$((FAILED_RUNS + 1))
                STATUS="FAILED"
                echo "STATUS: $STATUS (exit code: $EXIT_CODE)" | tee -a "$MAIN_LOG"
                echo "Failed after: ${ELAPSED}s ($ELAPSED_FORMATTED)" | tee -a "$MAIN_LOG"
                echo "Check logs:" | tee -a "$MAIN_LOG"
                echo "  stdout: $STDOUT_LOG" | tee -a "$MAIN_LOG"
                echo "  stderr: $STDERR_LOG" | tee -a "$MAIN_LOG"

                # Record timing to CSV
                echo "$TOTAL_RUNS,$genome,$n_beast,$max_class,$STATUS,$ELAPSED,$ELAPSED_MIN,$ELAPSED_FORMATTED,$START_TIME_FORMATTED,$END_TIME_FORMATTED" >> "$TIMING_CSV"
            fi
            echo "" | tee -a "$MAIN_LOG"
        done
    done

    echo "Finished processing genome: $genome" | tee -a "$MAIN_LOG"
    echo "" | tee -a "$MAIN_LOG"
done

# Final summary
echo "========================================" | tee -a "$MAIN_LOG"
echo "SUMMARY" | tee -a "$MAIN_LOG"
echo "========================================" | tee -a "$MAIN_LOG"
echo "Total runs: $TOTAL_RUNS" | tee -a "$MAIN_LOG"
echo "Successful: $SUCCESSFUL_RUNS" | tee -a "$MAIN_LOG"
echo "Failed: $FAILED_RUNS" | tee -a "$MAIN_LOG"
echo "Skipped: $SKIPPED_RUNS" | tee -a "$MAIN_LOG"
echo "Completed at: $(date)" | tee -a "$MAIN_LOG"
echo "" | tee -a "$MAIN_LOG"
echo "Output files:" | tee -a "$MAIN_LOG"
echo "  Main log: $MAIN_LOG" | tee -a "$MAIN_LOG"
echo "  Timing data (CSV): $TIMING_CSV" | tee -a "$MAIN_LOG"
echo "========================================" | tee -a "$MAIN_LOG"

if [ $FAILED_RUNS -gt 0 ]; then
    echo ""
    echo "Some runs failed. Check the log file for details: $MAIN_LOG"
    exit 1
fi