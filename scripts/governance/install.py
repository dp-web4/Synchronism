#!/usr/bin/env python3
"""
Installation Script for Synchronism Governance System

This script sets up the necessary environment for the self-governing repository management system.
"""

import os
import sys
import json
import subprocess
import argparse
from pathlib import Path

# Constants
# this script is run from Synchronism/scripts/governance so the paths are accordingly
REPO_PATH = os.getenv("REPO_PATH", os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

SCRIPTS_PATH = os.path.join(REPO_PATH, "scripts", "governance")
CONFIG_PATH = os.path.join(SCRIPTS_PATH, "config")

def create_directories():
    """Create necessary directories for the governance system."""
    print("Creating necessary directories...")
    
    # Create scripts directory
    os.makedirs(SCRIPTS_PATH, exist_ok=True)
    
    # Create config directory
    os.makedirs(CONFIG_PATH, exist_ok=True)
    
    print("Directories created successfully.")

def install_dependencies():
    """Install required Python dependencies."""
    print("Installing required dependencies...")
    
    # List of required packages
    requirements = [
        "pyyaml",
        "requests",
        "gitpython"
    ]
    
    # Install packages
    for package in requirements:
        print(f"Installing {package}...")
        try:
            subprocess.run([sys.executable, "-m", "pip", "install", package], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error installing {package}: {e}")
            return False
    
    print("Dependencies installed successfully.")
    return True

def setup_git_config():
    """Set up Git configuration for the governance system."""
    print("Setting up Git configuration...")
    
    try:
        # Check if Git is installed
        subprocess.run(["git", "--version"], check=True, capture_output=True)
        
        # Set Git configuration for GitHub Actions
        subprocess.run(["git", "config", "--local", "user.email", "actions@github.com"], check=True)
        subprocess.run(["git", "config", "--local", "user.name", "GitHub Actions"], check=True)
        
        print("Git configuration set up successfully.")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error setting up Git configuration: {e}")
        return False

def create_workflow_file():
    """Create GitHub Actions workflow file for daily updates."""
    print("Creating GitHub Actions workflow file...")
    
    workflows_dir = os.path.join(REPO_PATH, ".github", "workflows")
    os.makedirs(workflows_dir, exist_ok=True)
    
    workflow_path = os.path.join(workflows_dir, "synchronism_governance.yml")
    
    workflow_content = """name: Synchronism Governance System

on:
  schedule:
    - cron: '0 0 * * *'  # Run daily at midnight UTC
  workflow_dispatch:  # Allow manual triggering

jobs:
  governance_update:
    runs-on: ubuntu-latest
    
    steps:
    - name: Checkout repository
      uses: actions/checkout@v2
      
    - name: Set up Python
      uses: actions/setup-python@v2
      with:
        python-version: '3.10'
        
    - name: Install dependencies
      run: |
        python -m pip install --upgrade pip
        pip install pyyaml requests gitpython
        
    - name: Run governance system update
      run: |
        export REPO_PATH=$GITHUB_WORKSPACE
        python scripts/governance/main.py update
        
    - name: Commit and push changes
      run: |
        git config --local user.email "actions@github.com"
        git config --local user.name "GitHub Actions"
        git add .
        git diff --quiet && git diff --staged --quiet || git commit -m "Automated governance system update [$(date +'%Y-%m-%d')]"
        git push
"""
    
    with open(workflow_path, 'w') as f:
        f.write(workflow_content)
    
    print(f"Workflow file created at {workflow_path}")

def main():
    """Main function to run when script is executed directly."""
    parser = argparse.ArgumentParser(description="Install Synchronism Governance System")
    parser.add_argument("--repo-path", help="Path to the repository (default: current directory)")
    
    args = parser.parse_args()
    
    if args.repo_path:
        global REPO_PATH, SCRIPTS_PATH, CONFIG_PATH
        REPO_PATH = args.repo_path
        SCRIPTS_PATH = os.path.join(REPO_PATH, "scripts", "governance")
        CONFIG_PATH = os.path.join(SCRIPTS_PATH, "config")
    
    print("Installing Synchronism Governance System...")
    
    # Create directories
    create_directories()
    
    # Install dependencies
    if not install_dependencies():
        print("Error: Failed to install dependencies.")
        sys.exit(1)
    
    # Set up Git configuration
    if not setup_git_config():
        print("Warning: Failed to set up Git configuration.")
    
    # Create workflow file
    create_workflow_file()
    
    print("\nSynchronism Governance System installed successfully!")
    print("\nTo initialize the system, run:")
    print(f"  python {os.path.join(SCRIPTS_PATH, 'main.py')} init")
    print("\nTo run tests, run:")
    print(f"  python {os.path.join(SCRIPTS_PATH, 'test.py')}")
    print("\nTo run a daily update manually, run:")
    print(f"  python {os.path.join(SCRIPTS_PATH, 'main.py')} update")

if __name__ == "__main__":
    main()
