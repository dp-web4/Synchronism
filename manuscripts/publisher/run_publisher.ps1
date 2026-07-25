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

# Run Claude in WSL with the publisher prompt.
# Claude's stdout/stderr MUST be captured: without this redirect the agent's output goes to
# Task Scheduler's discarded stdout, so a failing run is indistinguishable from a successful
# one. Every log from 2026-07-04 to 2026-07-24 was exactly 270 bytes (header-only) — the
# publisher agent itself kept noting "cron 'Starting'-only" and covering the day manually.
#
# Do NOT use `claude -c`: it resumes the most recent conversation in this project store, and
# once that conversation outgrows the context limit every run dies at startup with
# "Prompt is too long" — self-perpetuating, because `-c` can never heal. This is the same
# root cause captured for the Archivist on 2026-07-15 (see private-context/archivist/
# run_archivist.ps1); that fix was applied to the sibling launcher only and never audited
# across the class. Each daily run must start a FRESH session; durable state lives in
# publisher/state/, publisher/log.md, and agent memory, not in conversation continuity.
wsl -d Ubuntu -u dp bash -c "cd /mnt/c/exe/projects/ai-agents/Synchronism/manuscripts && claude --dangerously-skip-permissions -p 'You are the Publisher. Read publisher/CLAUDE.md for your instructions. Execute the daily workflow: review SESSION_MAP for new complete arcs, evaluate publication candidates, update recommendations, track status changes, generate daily report. Log your activity to publisher/logs/.'" 2>&1 | Add-Content -Path $logFile

$claudeExit = $LASTEXITCODE

# Log completion
Add-Content -Path $logFile -Value "========================================="
Add-Content -Path $logFile -Value "Publisher Session Complete: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') (claude exit=$claudeExit)"
Add-Content -Path $logFile -Value "========================================="
