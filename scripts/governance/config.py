#!/usr/bin/env python3
"""
Configuration Script for Synchronism Governance System

This script creates initial configuration files for the governance system.
"""

import os
import sys
import json
import datetime
from pathlib import Path

# Constants
# this script is run from Synchronism/scripts/governance so the paths are accordingly
REPO_PATH = os.getenv("REPO_PATH", os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

SCRIPTS_PATH = os.path.join(REPO_PATH, "scripts", "governance")
CONFIG_PATH = os.path.join(SCRIPTS_PATH, "config")

# Fractal Scale Tags
FRACTAL_SCALES = {
    "quantum": "Quantum scale phenomena and models",
    "molecular": "Molecular scale interactions and structures",
    "biospheric": "Biological systems and ecosystems",
    "planetary": "Planetary scale processes and systems",
    "galactic": "Cosmic scale phenomena and models"
}

def create_config_directories():
    """Create necessary configuration directories."""
    print("Creating configuration directories...")
    
    # Create config directory
    os.makedirs(CONFIG_PATH, exist_ok=True)
    
    print("Configuration directories created successfully.")

def create_initial_config_files():
    """Create initial configuration files."""
    print("Creating initial configuration files...")
    
    # Create tokens.json
    tokens = {
        "contributors": {},
        "total_supply": {scale: 100 for scale in FRACTAL_SCALES},
        "last_distribution": datetime.datetime.now().isoformat()
    }
    
    with open(os.path.join(CONFIG_PATH, "tokens.json"), 'w') as f:
        json.dump(tokens, f, indent=2)
    
    # Create contributions.json
    contributions = {
        "contributions": []
    }
    
    with open(os.path.join(CONFIG_PATH, "contributions.json"), 'w') as f:
        json.dump(contributions, f, indent=2)
    
    # Create tensors.json
    tensors = {
        "contributors": {}
    }
    
    with open(os.path.join(CONFIG_PATH, "tensors.json"), 'w') as f:
        json.dump(tensors, f, indent=2)
    
    # Create reviews.json
    reviews = {
        "reviews": [],
        "review_assignments": {},
        "review_metrics": {
            "total_reviews": 0,
            "reviews_by_scale": {scale: 0 for scale in FRACTAL_SCALES},
            "reviews_by_type": {"human": 0, "ai": 0}
        }
    }
    
    with open(os.path.join(CONFIG_PATH, "reviews.json"), 'w') as f:
        json.dump(reviews, f, indent=2)
    
    # Create branches.json
    branches = {
        "branches": [
            {
                "name": "main",
                "parent": None,
                "scale": None,
                "type": None,
                "description": "Main branch of the Synchronism repository",
                "created_at": datetime.datetime.now().isoformat(),
                "status": "active",
                "children": []
            }
        ],
        "branch_metrics": {
            "total_branches": 1,
            "branches_by_scale": {scale: 0 for scale in FRACTAL_SCALES},
            "branches_by_type": {
                "exploration": 0,
                "refinement": 0,
                "integration": 0,
                "application": 0
            }
        }
    }
    
    with open(os.path.join(CONFIG_PATH, "branches.json"), 'w') as f:
        json.dump(branches, f, indent=2)
    
    # Create integration_log.json
    integration_log = {
        "integrations": []
    }
    
    with open(os.path.join(CONFIG_PATH, "integration_log.json"), 'w') as f:
        json.dump(integration_log, f, indent=2)
    
    # Create governance_log.json
    governance_log = {
        "events": [],
        "last_run": None,
        "system_status": "initialized",
        "version": "1.0.0"
    }
    
    with open(os.path.join(CONFIG_PATH, "governance_log.json"), 'w') as f:
        json.dump(governance_log, f, indent=2)
    
    print("Initial configuration files created successfully.")

def main():
    """Main function to run when script is executed directly."""
    print("Configuring Synchronism Governance System...")
    
    # Create directories
    create_config_directories()
    
    # Create initial config files
    create_initial_config_files()
    
    print("\nSynchronism Governance System configured successfully!")

if __name__ == "__main__":
    main()
