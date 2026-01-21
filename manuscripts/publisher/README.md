# Publisher Track Setup

Daily autonomous session that identifies publication-worthy research and tracks recommendations.

## Schedule

- **Time**: 02:30 UTC daily (1 hour after Archivist)
- **Machine**: CBP (Windows with WSL2)
- **Dependency**: Runs after Archivist updates SESSION_MAP

## Current Phase: 0 (Catalog and Recommend)

The Publisher currently:
- Catalogs existing publications
- Identifies publication-worthy session blocks
- Makes recommendations with rationale
- Tracks which recommendations are acted upon

Human performs actual publication decisions and execution.

## Windows Task Scheduler Setup

1. Open **Task Scheduler** (taskschd.msc)

2. Click **Create Basic Task**

3. **Name**: `Publisher Daily Session`
   **Description**: `Identifies publication-worthy research and tracks recommendations`

4. **Trigger**: Daily
   **Start**: 02:30:00
   **Recur every**: 1 day

5. **Action**: Start a program
   - **Program/script**: `powershell.exe`
   - **Arguments**: `-ExecutionPolicy Bypass -File "C:\exe\projects\ai-agents\Synchronism\manuscripts\publisher\run_publisher.ps1"`
   - **Start in**: `C:\exe\projects\ai-agents\Synchronism\manuscripts\publisher`

6. **Finish** and verify in Task Scheduler Library

## Manual Test

```powershell
# From PowerShell
powershell.exe -ExecutionPolicy Bypass -File "C:\exe\projects\ai-agents\Synchronism\manuscripts\publisher\run_publisher.ps1"
```

Or from WSL:
```bash
bash /mnt/c/exe/projects/ai-agents/Synchronism/manuscripts/publisher/run_publisher.sh
```

## Files

| File | Purpose |
|------|---------|
| `CLAUDE.md` | Publisher primer (evaluation criteria, workflow) |
| `run_publisher.ps1` | Windows launcher script |
| `run_publisher.sh` | WSL launcher script |
| `state/recommendations.json` | Publication recommendations |
| `state/published.json` | Catalog of published work |
| `logs/` | Session logs |
| `reports/` | Daily reports (created by Publisher) |

## Human Interaction

Review recommendations in `state/recommendations.json`:

```json
{
  "id": "REC-2026-001",
  "sessions": [285, 286, 287, 288, 289],
  "arc_name": "Quantum Computing Arc",
  "status": "recommended",
  "human_notes": ""  // Add your feedback here
}
```

Update status as you act on recommendations:
- `recommended` → `in_progress` → `published`
- Or `recommended` → `declined`

## Phase Roadmap

| Phase | Status | Capability |
|-------|--------|------------|
| 0 | **Current** | Catalog, recommend, track |
| 1 | Future | Draft preprints and abstracts |
| 2 | Future | Full preprint generation |

## Integration

- **Archivist** (01:30): Updates SESSION_MAP with latest sessions
- **Publisher** (02:30): Reviews SESSION_MAP, makes recommendations
- **Human**: Reviews recommendations, makes decisions, executes publication
