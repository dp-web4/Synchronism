# Proposal Lifecycle Management System

## Overview
Enhanced governance system with complete proposal lifecycle management, including withdrawal mechanisms, duplicate detection, archival, and automated maintenance.

## New Components

### 1. Proposal Cleanup (`proposal_cleanup.py`)
Manages cleanup of test proposals and archival of historical proposals.

**Features:**
- Automatic detection of test/duplicate proposals
- Intelligent archival decisions based on proposal value
- 90-day retention policy with permanent historical preservation option
- Comprehensive cleanup reports

**Usage:**
```bash
# Dry run (see what would be done)
python3 scripts/governance/proposal_cleanup.py --dry-run

# Execute cleanup
python3 scripts/governance/proposal_cleanup.py --execute
```

### 2. Enhanced Governance (`whitepaper_governance_enhanced.py`)
Extends the base governance system with lifecycle management features.

**New Proposal Statuses:**
- `withdrawn` - Author or arbiter withdrew the proposal
- `expired` - Automatically expired after 30 days of inactivity
- `duplicate` - Marked as duplicate of another proposal
- `superseded` - Replaced by a newer proposal

**Key Features:**
- **Withdrawal Mechanism**: Authors can withdraw proposals with reason tracking
- **Duplicate Detection**: 85% similarity threshold prevents redundant proposals
- **Auto-expiration**: Inactive proposals expire after 30 days
- **Supersession Tracking**: Links old proposals to their replacements

**API Examples:**
```python
from whitepaper_governance_enhanced import EnhancedWhitepaperGovernance, WithdrawalReason

governance = EnhancedWhitepaperGovernance()

# Withdraw a proposal
governance.withdraw_proposal(
    proposal_id="001",
    reason=WithdrawalReason.AUTHOR_REQUEST,
    withdrawer="Author Name",
    comment="Found a better approach"
)

# Check for duplicates before creating
duplicates = governance.check_duplicates(section, title, content)
if duplicates:
    print(f"Warning: Found {len(duplicates)} similar proposals")

# Expire old proposals
expired = governance.expire_old_proposals(dry_run=False)

# Supersede an old proposal
governance.supersede_proposal(
    old_id="001",
    new_id="002",
    superseder="Reviewer Name"
)
```

### 3. Governance Maintenance (`governance_maintenance.py`)
Automated maintenance system for ongoing governance health.

**Maintenance Tasks:**
1. Expire proposals inactive for 30+ days
2. Detect and report duplicate proposals
3. Archive completed proposals after 7 days
4. Generate comprehensive maintenance reports

**Usage:**
```bash
# Dry run maintenance
python3 scripts/governance/governance_maintenance.py --dry-run

# Execute maintenance
python3 scripts/governance/governance_maintenance.py --execute

# Execute and save report
python3 scripts/governance/governance_maintenance.py --execute --save-report
```

## Archive System

### Structure
```json
{
  "archived_proposals": [
    {
      "id": "001",
      "title": "Proposal Title",
      "status": "rejected",
      "archive_reason": "rejected",
      "archive_date": "2025-08-21T10:00:00",
      "preserve_until": "2025-11-21T10:00:00",  // or null for permanent
      "historical_value": true
    }
  ],
  "archive_policy": {
    "retention_days": 90,
    "last_purge": "2025-08-21",
    "total_archived": 15,
    "total_purged": 0
  }
}
```

### Archive Decisions

**Always Archive (with historical_value=true):**
- Accepted/implemented proposals
- Proposals with meaningful reviews
- Superseded proposals showing evolution

**Archive Temporarily (90 days):**
- Withdrawn proposals after discussion
- Expired proposals with some activity
- Rejected proposals without extensive review

**Never Archive (Delete):**
- Test proposals
- Duplicates without reviews
- Malformed/spam proposals
- Immediately withdrawn mistakes

## Integration with Daily Governance

The enhanced system integrates seamlessly with the existing daily governance cycle:

1. **Pre-cycle Maintenance** (optional):
   ```python
   maintenance = GovernanceMaintenance()
   maintenance.run_maintenance(dry_run=False)
   ```

2. **Duplicate Check During Creation**:
   - Automatically warns about similar proposals
   - Prevents accidental duplicates

3. **Lifecycle Management**:
   - Authors can withdraw within grace period
   - Old proposals auto-expire
   - Completed proposals auto-archive

4. **Post-cycle Cleanup** (weekly recommended):
   ```bash
   python3 scripts/governance/governance_maintenance.py --execute
   ```

## Configuration

### Default Settings
```python
# In whitepaper_governance_enhanced.py
expiration_days = 30         # Days before auto-expiration
similarity_threshold = 0.85   # Duplicate detection threshold

# In proposal_cleanup.py
retention_days = 90          # Archive retention period
```

### Customization
Settings can be modified in the respective Python files or passed as parameters to the class constructors.

## Workflow Examples

### Handle Duplicate Test Proposals
```bash
# 1. Initial cleanup of test proposals
python3 scripts/governance/proposal_cleanup.py --execute

# 2. Check results
cat scripts/governance/config/whitepaper_proposals.json
cat scripts/governance/config/whitepaper_proposals_archive.json
```

### Author Withdraws Proposal
```python
governance = EnhancedWhitepaperGovernance()
result = governance.withdraw_proposal(
    proposal_id="002",
    reason=WithdrawalReason.AUTHOR_REQUEST,
    withdrawer="Claude-3.5",
    comment="Realized this conflicts with section 5"
)
```

### Weekly Maintenance
```bash
# Run every Sunday via cron or GitHub Actions
python3 scripts/governance/governance_maintenance.py --execute --save-report
```

## Benefits

1. **Clean Active List**: Only relevant proposals remain active
2. **Historical Preservation**: Important discussions preserved in archive
3. **Reduced Redundancy**: Duplicate detection prevents wasted effort
4. **Automated Hygiene**: System self-maintains without manual intervention
5. **Clear Lifecycle**: Proposals have defined states and transitions
6. **Audit Trail**: All actions tracked with timestamps and reasons

## Migration from Old System

The system is backward compatible. Existing proposals are preserved and can be managed with the new lifecycle features. The initial cleanup (`proposal_cleanup.py --execute`) migrated test proposals to the archive while preserving legitimate proposals.

## Future Enhancements

Potential improvements for consideration:
- Configurable expiration periods per proposal type
- Batch operations UI for easier management
- Metrics dashboard for proposal lifecycle analytics
- Integration with notification system for withdrawals/expirations
- Automatic duplicate merging suggestions

---

*Last Updated: August 21, 2025*
*System Version: 2.0.0*