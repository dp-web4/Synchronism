#!/bin/bash
# Persistent simulation runner using nohup with proper process management
# Creates a managed process that survives session termination

set -e

USAGE="Usage: $0 <job_name> <python_script> [args...]
Example: $0 su2_physics synchronism_session30_su2_lattice_3p1d.py"

if [ "$#" -lt 2 ]; then
    echo "$USAGE"
    exit 1
fi

JOB_NAME="$1"
SCRIPT="$2"
shift 2
ARGS="$@"

# Directories
SCRIPT_DIR="/mnt/c/exe/projects/ai-agents/synchronism/simulations"
PID_DIR="$SCRIPT_DIR/.pids"
LOG_DIR="$SCRIPT_DIR"

mkdir -p "$PID_DIR"

SCRIPT_PATH="$SCRIPT_DIR/$SCRIPT"
PID_FILE="$PID_DIR/${JOB_NAME}.pid"
LOG_FILE="$LOG_DIR/${JOB_NAME}.log"
STATUS_FILE="$PID_DIR/${JOB_NAME}.status"

# Check if script exists
if [ ! -f "$SCRIPT_PATH" ]; then
    echo "Error: Script not found: $SCRIPT_PATH"
    exit 1
fi

# Check if job already running
if [ -f "$PID_FILE" ]; then
    OLD_PID=$(cat "$PID_FILE")
    if ps -p "$OLD_PID" > /dev/null 2>&1; then
        echo "Error: Job '$JOB_NAME' already running (PID: $OLD_PID)"
        echo "Log: tail -f $LOG_FILE"
        echo "Kill: kill $OLD_PID"
        exit 1
    else
        echo "Removing stale PID file..."
        rm -f "$PID_FILE"
    fi
fi

# Launch job
echo "Launching persistent job: $JOB_NAME"
echo "Script: $SCRIPT $ARGS"
echo "Log: $LOG_FILE"
echo ""

cd "$SCRIPT_DIR"

# Write status
echo "STARTED $(date)" > "$STATUS_FILE"

# Launch with nohup, redirect all output, run in background
nohup python3 "$SCRIPT" $ARGS > "$LOG_FILE" 2>&1 &
JOB_PID=$!

# Save PID
echo "$JOB_PID" > "$PID_FILE"

# Wait a moment to check if it started
sleep 3

if ps -p "$JOB_PID" > /dev/null 2>&1; then
    echo "✓ Job started successfully"
    echo "  PID: $JOB_PID"
    echo "  Log: $LOG_FILE"
    echo ""
    echo "Monitor:"
    echo "  tail -f $LOG_FILE"
    echo ""
    echo "Check status:"
    echo "  ps -p $JOB_PID"
    echo ""
    echo "Kill job:"
    echo "  kill $JOB_PID"
    echo "  # or"
    echo "  $0 stop $JOB_NAME"

    # Update status
    echo "RUNNING $(date) PID=$JOB_PID" >> "$STATUS_FILE"
else
    echo "✗ Job failed to start"
    rm -f "$PID_FILE"
    echo "FAILED $(date)" >> "$STATUS_FILE"
    exit 1
fi
