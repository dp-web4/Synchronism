#!/bin/bash
# Screen-based simulation runner for long-running Synchronism simulations
# Solves infrastructure issue: background processes terminating when sessions end

set -e

USAGE="Usage: $0 <screen_name> <python_script> [args...]
Example: $0 su2_physics synchronism_session30_su2_lattice_3p1d.py"

if [ "$#" -lt 2 ]; then
    echo "$USAGE"
    exit 1
fi

SCREEN_NAME="$1"
SCRIPT="$2"
shift 2
ARGS="$@"

# Full script path
SCRIPT_DIR="/mnt/c/exe/projects/ai-agents/synchronism/simulations"
SCRIPT_PATH="$SCRIPT_DIR/$SCRIPT"

# Check if script exists
if [ ! -f "$SCRIPT_PATH" ]; then
    echo "Error: Script not found: $SCRIPT_PATH"
    exit 1
fi

# Check if screen session already exists
if screen -list | grep -q "\.$SCREEN_NAME\s"; then
    echo "Error: Screen session '$SCREEN_NAME' already exists"
    echo "To attach: screen -r $SCREEN_NAME"
    echo "To kill: screen -S $SCREEN_NAME -X quit"
    exit 1
fi

# Create detached screen session
echo "Creating screen session: $SCREEN_NAME"
echo "Running: python3 $SCRIPT $ARGS"
echo ""

cd "$SCRIPT_DIR"
screen -dmS "$SCREEN_NAME" bash -c "python3 $SCRIPT $ARGS 2>&1 | tee ${SCREEN_NAME}.log; echo 'Simulation complete. Press Enter to close.'; read"

# Wait a moment for screen to start
sleep 2

# Check if session is running
if screen -list | grep -q "\.$SCREEN_NAME\s"; then
    echo "✓ Screen session '$SCREEN_NAME' created successfully"
    echo ""
    echo "Monitor progress:"
    echo "  tail -f $SCRIPT_DIR/${SCREEN_NAME}.log"
    echo ""
    echo "Attach to session:"
    echo "  screen -r $SCREEN_NAME"
    echo ""
    echo "List all sessions:"
    echo "  screen -ls"
    echo ""
    echo "Kill session:"
    echo "  screen -S $SCREEN_NAME -X quit"
else
    echo "✗ Failed to create screen session"
    exit 1
fi
