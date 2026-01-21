# PowerShell script to start Publisher session
# Run from Windows Task Scheduler at 02:30 daily
#
# Task Scheduler Setup:
#   1. Open Task Scheduler
#   2. Create Basic Task: "Publisher Daily Session"
#   3. Trigger: Daily at 02:30
#   4. Action: Start a program
#   5. Program: powershell.exe
#   6. Arguments: -ExecutionPolicy Bypass -File "C:\exe\projects\ai-agents\Synchronism\manuscripts\publisher\run_publisher.ps1"
#   7. Start in: C:\exe\projects\ai-agents\Synchronism\manuscripts

$logFile = "C:\exe\projects\ai-agents\Synchronism\manuscripts\publisher\logs\publisher-$(Get-Date -Format 'yyyy-MM-dd').log"

# Log start
Add-Content -Path $logFile -Value "========================================="
Add-Content -Path $logFile -Value "Publisher Session Starting: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')"
Add-Content -Path $logFile -Value "========================================="

# Run Claude in WSL with the publisher prompt
wsl -d Ubuntu -u dp bash -c "cd /mnt/c/exe/projects/ai-agents/Synchronism/manuscripts && claude -c --dangerously-skip-permissions -p 'You are the Publisher. Read publisher/CLAUDE.md for your instructions. Execute the daily workflow: review SESSION_MAP for new complete arcs, evaluate publication candidates, update recommendations, track status changes, generate daily report. Log your activity to publisher/logs/.'"

# Log completion
Add-Content -Path $logFile -Value "========================================="
Add-Content -Path $logFile -Value "Publisher Session Complete: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')"
Add-Content -Path $logFile -Value "========================================="
