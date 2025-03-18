#!/usr/bin/env python3
"""
Fractal Branch Management Script for Synchronism Governance System

This script manages the fractal branch structure of the Synchronism repository,
allowing focused explorations of different aspects of reality within the model.
"""

import os
import sys
import json
import datetime
import subprocess
import re
from pathlib import Path

# Constants
REPO_PATH = os.getenv("REPO_PATH", os.getcwd())
CONFIG_PATH = os.path.join(REPO_PATH, "scripts", "governance", "config")
BRANCHES_PATH = os.path.join(CONFIG_PATH, "branches.json")

# Fractal Scale Tags
FRACTAL_SCALES = {
    "quantum": "Quantum scale phenomena and models",
    "molecular": "Molecular scale interactions and structures",
    "biospheric": "Biological systems and ecosystems",
    "planetary": "Planetary scale processes and systems",
    "galactic": "Cosmic scale phenomena and models"
}

# Branch Types
BRANCH_TYPES = {
    "exploration": "Exploratory branches for new concepts",
    "refinement": "Refinement branches for existing concepts",
    "integration": "Integration branches for merging concepts",
    "application": "Application branches for practical implementations"
}

class FractalBranchSystem:
    """Manages the fractal branch structure of the Synchronism repository."""
    
    def __init__(self):
        """Initialize the fractal branch system."""
        self.branches = self._load_branches()
        
    def _load_branches(self):
        """Load branch data from file or create if not exists."""
        os.makedirs(CONFIG_PATH, exist_ok=True)
        
        if os.path.exists(BRANCHES_PATH):
            with open(BRANCHES_PATH, 'r') as f:
                return json.load(f)
        else:
            # Initialize with empty branch data and main branch
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
                    "branches_by_type": {btype: 0 for btype in BRANCH_TYPES}
                }
            }
            self._save_branches(branches)
            return branches
    
    def _save_branches(self, branches=None):
        """Save branch data to file."""
        if branches is None:
            branches = self.branches
            
        with open(BRANCHES_PATH, 'w') as f:
            json.dump(branches, f, indent=2)
    
    def get_branch(self, branch_name):
        """Get branch data by name."""
        for branch in self.branches["branches"]:
            if branch["name"] == branch_name:
                return branch
        return None
    
    def create_branch(self, name, parent_name, scale, branch_type, description):
        """
        Create a new branch in the fractal structure.
        
        Args:
            name: Name of the new branch
            parent_name: Name of the parent branch
            scale: Scale tag for the branch
            branch_type: Type of branch
            description: Description of the branch purpose
            
        Returns:
            dict: Result of the branch creation process
        """
        # Validate inputs
        if scale not in FRACTAL_SCALES:
            return {
                "status": "error",
                "message": f"Invalid scale: {scale}. Must be one of {list(FRACTAL_SCALES.keys())}"
            }
        
        if branch_type not in BRANCH_TYPES:
            return {
                "status": "error",
                "message": f"Invalid branch type: {branch_type}. Must be one of {list(BRANCH_TYPES.keys())}"
            }
        
        # Check if branch already exists
        if self.get_branch(name):
            return {
                "status": "error",
                "message": f"Branch already exists: {name}"
            }
        
        # Get parent branch
        parent_branch = self.get_branch(parent_name)
        if not parent_branch:
            return {
                "status": "error",
                "message": f"Parent branch not found: {parent_name}"
            }
        
        # Create branch in Git
        try:
            # Make sure we're on the parent branch
            subprocess.run(["git", "checkout", parent_name], check=True)
            
            # Create new branch
            subprocess.run(["git", "checkout", "-b", name], check=True)
            
            # Create branch metadata file
            metadata = {
                "name": name,
                "parent": parent_name,
                "scale": scale,
                "type": branch_type,
                "description": description,
                "created_at": datetime.datetime.now().isoformat()
            }
            
            metadata_path = os.path.join(REPO_PATH, f"branch-{name}.json")
            with open(metadata_path, 'w') as f:
                json.dump(metadata, f, indent=2)
            
            # Commit the metadata file
            subprocess.run(["git", "add", metadata_path], check=True)
            subprocess.run(["git", "commit", "-m", f"Create branch metadata for {name}"], check=True)
            
            # Push the branch
            subprocess.run(["git", "push", "--set-upstream", "origin", name], check=True)
            
            # Return to main branch
            subprocess.run(["git", "checkout", "main"], check=True)
            
        except subprocess.CalledProcessError as e:
            return {
                "status": "error",
                "message": f"Git operation failed: {str(e)}",
                "details": e.stderr if hasattr(e, 'stderr') else None
            }
        
        # Create branch record
        branch = {
            "name": name,
            "parent": parent_name,
            "scale": scale,
            "type": branch_type,
            "description": description,
            "created_at": datetime.datetime.now().isoformat(),
            "status": "active",
            "children": []
        }
        
        # Add to branches list
        self.branches["branches"].append(branch)
        
        # Add as child to parent
        parent_branch["children"].append(name)
        
        # Update metrics
        self.branches["branch_metrics"]["total_branches"] += 1
        self.branches["branch_metrics"]["branches_by_scale"][scale] += 1
        self.branches["branch_metrics"]["branches_by_type"][branch_type] += 1
        
        # Save changes
        self._save_branches()
        
        return {
            "status": "success",
            "message": f"Branch {name} created successfully",
            "branch": branch
        }
    
    def merge_branch(self, source_name, target_name):
        """
        Merge a source branch into a target branch.
        
        Args:
            source_name: Name of the source branch
            target_name: Name of the target branch
            
        Returns:
            dict: Result of the merge process
        """
        # Get branch data
        source_branch = self.get_branch(source_name)
        target_branch = self.get_branch(target_name)
        
        if not source_branch:
            return {
                "status": "error",
                "message": f"Source branch not found: {source_name}"
            }
        
        if not target_branch:
            return {
                "status": "error",
                "message": f"Target branch not found: {target_name}"
            }
        
        # Perform Git merge
        try:
            # Make sure we're on the target branch
            subprocess.run(["git", "checkout", target_name], check=True)
            
            # Merge source branch
            subprocess.run(["git", "merge", source_name], check=True)
            
            # Push the changes
            subprocess.run(["git", "push"], check=True)
            
            # Return to main branch
            subprocess.run(["git", "checkout", "main"], check=True)
            
        except subprocess.CalledProcessError as e:
            return {
                "status": "error",
                "message": f"Git operation failed: {str(e)}",
                "details": e.stderr if hasattr(e, 'stderr') else None
            }
        
        # Update branch status
        source_branch["status"] = "merged"
        source_branch["merged_into"] = target_name
        source_branch["merged_at"] = datetime.datetime.now().isoformat()
        
        # Save changes
        self._save_branches()
        
        return {
            "status": "success",
            "message": f"Branch {source_name} merged into {target_name} successfully"
        }
    
    def archive_branch(self, branch_name):
        """
        Archive a branch that is no longer active.
        
        Args:
            branch_name: Name of the branch to archive
            
        Returns:
            dict: Result of the archive process
        """
        # Get branch data
        branch = self.get_branch(branch_name)
        
        if not branch:
            return {
                "status": "error",
                "message": f"Branch not found: {branch_name}"
            }
        
        # Check if branch has unmerged children
        active_children = []
        for child_name in branch["children"]:
            child = self.get_branch(child_name)
            if child and child["status"] == "active":
                active_children.append(child_name)
        
        if active_children:
            return {
                "status": "error",
                "message": f"Branch {branch_name} has active children: {active_children}"
            }
        
        # Update branch status
        branch["status"] = "archived"
        branch["archived_at"] = datetime.datetime.now().isoformat()
        
        # Save changes
        self._save_branches()
        
        return {
            "status": "success",
            "message": f"Branch {branch_name} archived successfully"
        }
    
    def initialize_fractal_structure(self):
        """
        Initialize the fractal branch structure with scale-specific branches.
        
        Returns:
            dict: Result of the initialization process
        """
        results = {
            "status": "success",
            "message": "Fractal branch structure initialized successfully",
            "created_branches": []
        }
        
        # Create scale branches
        for scale in FRACTAL_SCALES:
            branch_name = f"scale-{scale}"
            description = f"Branch for {FRACTAL_SCALES[scale]}"
            
            result = self.create_branch(
                name=branch_name,
                parent_name="main",
                scale=scale,
                branch_type="exploration",
                description=description
            )
            
            if result["status"] == "success":
                results["created_branches"].append(branch_name)
            else:
                # If branch already exists, continue
                if "already exists" in result["message"]:
                    continue
                
                # Otherwise, return error
                return {
                    "status": "error",
                    "message": f"Failed to create branch {branch_name}: {result['message']}"
                }
        
        return results
    
    def generate_branch_report(self):
        """Generate a report of the fractal branch structure."""
        report = {
            "timestamp": datetime.datetime.now().isoformat(),
            "total_branches": self.branches["branch_metrics"]["total_branches"],
            "branches_by_scale": self.branches["branch_metrics"]["branches_by_scale"],
            "branches_by_type": self.branches["branch_metrics"]["branches_by_type"],
            "active_branches": 0,
            "merged_branches": 0,
            "archived_branches": 0,
            "branch_tree": self._generate_branch_tree("main")
        }
        
        # Count branch statuses
        for branch in self.branches["branches"]:
            status = branch.get("status", "active")
            if status == "active":
                report["active_branches"] += 1
            elif status == "merged":
                report["merged_branches"] += 1
            elif status == "archived":
                report["archived_branches"] += 1
        
        # Save report
        report_path = os.path.join(CONFIG_PATH, f"branch_report_{datetime.datetime.now().strftime('%Y%m%d')}.json")
        with open(report_path, 'w') as f:
            json.dump(report, f, indent=2)
        
        return {
            "status": "success",
            "message": "Branch report generated successfully",
            "report_path": report_path,
            "report": report
        }
    
    def _generate_branch_tree(self, branch_name, depth=0):
        """
        Generate a tree representation of the branch structure.
        
        Args:
            branch_name: Name of the branch to start from
            depth: Current depth in the tree
            
        Returns:
            dict: Tree representation of the branch structure
        """
        branch = self.get_branch(branch_name)
        if not branch:
            return None
        
        tree = {
            "name": branch["name"],
            "scale": branch["scale"],
            "type": branch["type"],
            "status": branch["status"],
            "depth": depth,
            "children": []
        }
        
        for child_name in branch["children"]:
            child_tree = self._generate_branch_tree(child_name, depth + 1)
            if child_tree:
                tree["children"].append(child_tree)
        
        return tree


def main():
    """Main function to run when script is executed directly."""
    print("Synchronism Fractal Branch System")
    print("================================")
    
    # Create necessary directories
    os.makedirs(CONFIG_PATH, exist_ok=True)
    
    branch_system = FractalBranchSystem()
    
    # Initialize fractal structure if not already done
    if branch_system.branches["branch_metrics"]["total_branches"] <= 1:
        print("Initializing fractal branch structure...")
        init_result = branch_system.initialize_fractal_structure()
        print(f"Initialization result: {init_result['status']} - {init_result['message']}")
        if "created_branches" in init_result:
            print(f"Created branches: {init_result['created_branches']}")
    else:
        print(f"Fractal structure already initialized with {branch_system.branches['branch_metrics']['total_branches']} branches")
    
    # Generate branch report
    report_result = branch_system.generate_branch_report()
    print(f"Report result: {report_result['status']} - {report_result['message']}")
    print(f"Report saved to: {report_result['report_path']}")


if __name__ == "__main__":
    main()
