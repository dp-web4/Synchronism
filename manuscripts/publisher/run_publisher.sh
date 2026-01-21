#!/bin/bash
# Publisher daily session launcher
# Called by Windows Task Scheduler via PowerShell wrapper
#
# Direct invocation (for testing):
#   bash /mnt/c/exe/projects/ai-agents/Synchronism/manuscripts/publisher/run_publisher.sh

LOG_DIR="/mnt/c/exe/projects/ai-agents/Synchronism/manuscripts/publisher/logs"
LOG_FILE="$LOG_DIR/publisher-$(date +%Y-%m-%d).log"

echo "=========================================" >> "$LOG_FILE"
echo "Publisher Session Starting: $(date)" >> "$LOG_FILE"
echo "=========================================" >> "$LOG_FILE"

cd /mnt/c/exe/projects/ai-agents/Synchronism/manuscripts

# Run Claude with publisher instructions
claude -c --dangerously-skip-permissions -p "You are the Publisher. Read publisher/CLAUDE.md for your instructions. Execute the daily workflow: review SESSION_MAP for new complete arcs, evaluate publication candidates, update recommendations, track status changes, generate daily report. Log your activity to publisher/logs/." >> "$LOG_FILE" 2>&1

echo "=========================================" >> "$LOG_FILE"
echo "Publisher Session Complete: $(date)" >> "$LOG_FILE"
echo "=========================================" >> "$LOG_FILE"
