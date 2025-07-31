#!/bin/bash

set -e

if [ "$#" -lt 4 ]; then
    echo "Usage: $0 <input_file> <timeout> [component_params...]"
    exit 1
fi

EXE_PATH=$1
INPUT_FILE=$2
TIMEOUT=$3
GRID_SIZE=$4
ERASE_ALL_DUPLICATE=$5

# Use /usr/bin/time for memory usage (maximum resident set size in KB)
TMP_LOG=$(mktemp)

# Run the benchmarked command
/usr/bin/time -f "TIME:%e\nMEM:%M" timeout "$TIMEOUT"s "$EXE_PATH/performance_snap_polygon_soup" "$INPUT_FILE" "$GRID_SIZE" "$ERASE_ALL_DUPLICATE" 2> "$TMP_LOG"

# Parse time and memory
RUNNING_TIME=$(grep "TIME" "$TMP_LOG" | cut -d':' -f2)
MEMORY=$(grep "MEM" "$TMP_LOG" | cut -d':' -f2)

rm -f "$TMP_LOG"

# Output JSON
echo "{\"seconds\": \"$RUNNING_TIME\", \"memory_peaks\": \"$MEMORY\"}"
