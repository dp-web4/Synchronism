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
```

### Module Usage

- **Fractal Branch System**: Primary owner, creates and updates branch structure
- **Contribution System**: Uses branch information to assign contributions to scales
- **Integration System**: Uses branch information for integration operations

## contributions.json

### Purpose
Tracks all contributions to the repository, including their status, scale, and associated metadata.

### Structure
```json
{
  "contributions": [
    {
      "id": "contribution_id",
      "contributor_id": "contributor_id",
      "timestamp": "submission_timestamp",
      "content": "contribution_content",
      "file_path": "affected_file_path",
      "scale": "fractal_scale",
      "status": "contribution_status",
      "reviews": [],
      "validations": [],
      "value_metrics": {}
    }
  ]
}
```

### Module Usage

- **Contribution System**: Creates and updates contribution records
- **Review System**: Updates contribution status based on reviews
- **Validation System**: Updates contribution with validation results
- **Integration System**: Uses contribution data for integration
- **Token System**: Uses contribution data for token exchanges

## reviews.json

### Purpose
Records all reviews and review assignments, tracking review metrics and recommendations.

### Structure
```json
{
  "reviews": [
    {
      "id": "review_id",
      "contribution_id": "contribution_id",
      "reviewer_id": "reviewer_id",
      "timestamp": "review_timestamp",
      "recommendation": "review_recommendation",
      "comments": "review_comments",
      "metrics": {
        "value_score": 0.8,
        "veracity_score": 0.7,
        "validity_score": 0.9
      }
    }
  ],
  "review_assignments": {
    "contribution_id": ["reviewer_id_1", "reviewer_id_2"]
  },
  "review_metrics": {
    "total_reviews": 0,
    "reviews_by_scale": {
      "quantum": 0,
      "molecular": 0,
      "biospheric": 0,
      "planetary": 0,
      "galactic": 0
    },
    "reviews_by_type": {
      "human": 0,
      "ai": 0
    }
  }
}
```

### Module Usage

- **Review System**: Primary owner, creates and updates review records
- **Validation System**: Uses review data to inform validation
- **Token System**: Uses review data for reviewer token management

## tensors.json

### Purpose
Maintains contributor trust tensors using the T3/V3 framework, tracking talent, training, and temperament scores.

### Structure
```json
{
  "contributors": {
    "contributor_id": {
      "trust_tensor": {
        "talent": {
          "quantum": 0.5,
          "molecular": 0.5,
          "biospheric": 0.5,
          "planetary": 0.5,
          "galactic": 0.5
        },
        "training": {
          "quantum": 0.5,
          "molecular": 0.5,
          "biospheric": 0.5,
          "planetary": 0.5,
          "galactic": 0.5
        },
        "temperament": {
          "quality": 0.5,
          "accuracy": 0.5,
          "collaboration": 0.5
        }
      },
      "value_tensor": {
        "value": 0.5,
        "veracity": 0.5,
        "validity": 0.5
      },
      "contribution_count": 0,
      "last_updated": "timestamp"
    }
  }
}
```

### Module Usage

- **Validation System**: Primary owner, updates trust and value tensors
- **Review System**: Uses trust tensors for reviewer selection
- **Token System**: Uses trust tensors for token distribution

## tokens.json

### Purpose
Tracks token distribution and balances for all contributors, implementing the ATP/ADP-like token system.

### Structure
```json
{
  "contributors": {
    "contributor_id": {
      "charged": {
        "quantum": 10,
        "molecular": 10,
        "biospheric": 10,
        "planetary": 10,
        "galactic": 10
      },
      "discharged": {
        "quantum": 0,
        "molecular": 0,
        "biospheric": 0,
        "planetary": 0,
        "galactic": 0
      },
      "last_distribution": "timestamp"
    }
  },
  "total_supply": {
    "quantum": 100,
    "molecular": 100,
    "biospheric": 100,
    "planetary": 100,
    "galactic": 100
  },
  "last_distribution": "timestamp"
}
```

### Module Usage

- **Token System**: Primary owner, manages token distribution and exchanges
- **Contribution System**: Discharges tokens for contributions
- **Review System**: Manages tokens for review activities
- **Validation System**: Triggers token recharging based on validation

## Report Files

### Purpose
Generated periodically to document system status and metrics, providing snapshots of system activity.

### Types
- **branch_report_YYYYMMDD.json**: Documents branch structure and metrics
- **review_report_YYYYMMDD.json**: Documents review activity and metrics
- **token_report_YYYYMMDD.json**: Documents token distribution and usage
- **governance_log.json**: Tracks overall governance system activity

### Structure (Example: branch_report)
```json
{
  "timestamp": "report_generation_timestamp",
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
  },
  "active_branches": 6,
  "merged_branches": 0,
  "archived_branches": 0,
  "branch_tree": {}
}
```

### Module Usage

- **Fractal Branch System**: Generates branch reports
- **Review System**: Generates review reports
- **Token System**: Generates token reports
- **Main System**: Uses reports for system monitoring

## Integration Log (integration_log.json)

### Purpose
Maintains a record of all integrations performed by the system, tracking which contributions have been integrated.

### Structure
```json
{
  "integrations": [
    {
      "contribution_id": "contribution_id",
      "timestamp": "integration_timestamp",
      "scale": "fractal_scale",
      "value_score": 1.2
    }
  ]
}
```

### Module Usage
- **Integration System**: Primary owner, records integration activities
- **Main System**: Uses integration log for system monitoring

## Conclusion
The configuration files form the backbone of the Synchronism Governance System, storing all the data needed for the system to function. Each file has a specific purpose and is used by multiple modules, creating a web of interactions that enables the system to operate as a cohesive whole. The JSON format provides a flexible, human-readable structure that can be easily extended as the system evolves.