# Configuration Files

## Introduction

The Synchronism Governance System uses a set of JSON configuration files to maintain state, track metrics, and store system data. This document provides detailed documentation for each configuration file, explaining its purpose, structure, and which modules interact with it.

## Configuration Directory Structure

All configuration files are stored in the `scripts/governance/config` directory. The system uses the following files:

- **branches.json**: Stores the fractal branch structure
- **contributions.json**: Tracks all contributions and their status
- **reviews.json**: Records reviews and review assignments
- **tensors.json**: Maintains contributor trust tensors
- **tokens.json**: Tracks token distribution and balances
- **Report files**: Generated periodically to document system status

## branches.json

### Purpose
Stores the complete fractal branch structure of the repository, including branch relationships, scales, and metrics.

### Structure
```json
{
  "branches": [
    {
      "name": "branch_name",
      "parent": "parent_branch_name",
      "scale": "fractal_scale",
      "type": "branch_type",
      "description": "branch_description",
      "created_at": "timestamp",
      "status": "branch_status",
      "children": ["child_branch_1", "child_branch_2"]
    }
  ],
  "branch_metrics": {
    "total_branches": 6,
    "branches_by_scale": {
      "quantum": 1,
      "molecular": 1,
      "biospheric": 1,
      "planetary": 1,
      "galactic": 1
    },
    "branches_by_type": {
      "exploration": 5,
      "refinement": 0,
      "integration": 0,
      "application": 0
    }
  }
}
