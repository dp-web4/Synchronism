#!/usr/bin/env python3
"""
Integration Script for Synchronism Governance System

This script handles the integration of validated contributions into the Synchronism repository.
It implements the final stage of the contribution workflow, applying approved changes to the repository.
"""

import os
import sys
import json
import datetime
import subprocess
import re
from pathlib import Path

# Constants
# this script is run from Synchronism/scripts/governance so the paths are accordingly
REPO_PATH = os.getenv("REPO_PATH", os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

CONFIG_PATH = os.path.join(REPO_PATH, "scripts", "governance", "config")
CONTRIBUTIONS_PATH = os.path.join(CONFIG_PATH, "contributions.json")
INTEGRATION_LOG_PATH = os.path.join(CONFIG_PATH, "integration_log.json")

# Fractal Scale Tags
FRACTAL_SCALES = {
    "quantum": "Quantum scale phenomena and models",
    "molecular": "Molecular scale interactions and structures",
    "biospheric": "Biological systems and ecosystems",
    "planetary": "Planetary scale processes and systems",
    "galactic": "Cosmic scale phenomena and models"
}

class IntegrationSystem:
    """Integrates validated contributions into the Synchronism repository."""
    
    def __init__(self):
        """Initialize the integration system."""
        self.contributions = self._load_contributions()
        self.integration_log = self._load_integration_log()
        
    def _load_contributions(self):
        """Load contribution data from file."""
        if os.path.exists(CONTRIBUTIONS_PATH):
            with open(CONTRIBUTIONS_PATH, 'r') as f:
                return json.load(f)
        else:
            return {"contributions": []}
    
    def _save_contributions(self, contributions=None):
        """Save contribution data to file."""
        if contributions is None:
            contributions = self.contributions
            
        with open(CONTRIBUTIONS_PATH, 'w') as f:
            json.dump(contributions, f, indent=2)
    
    def _load_integration_log(self):
        """Load integration log from file or create if not exists."""
        os.makedirs(CONFIG_PATH, exist_ok=True)
        
        if os.path.exists(INTEGRATION_LOG_PATH):
            with open(INTEGRATION_LOG_PATH, 'r') as f:
                return json.load(f)
        else:
            # Initialize with empty log
            log = {
                "integrations": []
            }
            self._save_integration_log(log)
            return log
    
    def _save_integration_log(self, log=None):
        """Save integration log to file."""
        if log is None:
            log = self.integration_log
            
        with open(INTEGRATION_LOG_PATH, 'w') as f:
            json.dump(log, f, indent=2)
    
    def get_certified_contributions(self):
        """Get all certified contributions that are ready for integration."""
        return [c for c in self.contributions["contributions"] 
                if c["status"] == "certified" and c.get("value_certified", False)]
    
    def integrate_contribution(self, contribution_id):
        """
        Integrate a certified contribution into the repository.
        
        Args:
            contribution_id: Unique identifier for the contribution
            
        Returns:
            dict: Result of the integration process
        """
        # Find the contribution
        contribution = None
        for c in self.contributions["contributions"]:
            if c["id"] == contribution_id:
                contribution = c
                break
                
        if contribution is None:
            return {
                "status": "error",
                "message": f"Contribution not found: {contribution_id}"
            }
        
        # Check if contribution is certified
        if contribution["status"] != "certified" or not contribution.get("value_certified", False):
            return {
                "status": "error",
                "message": f"Contribution {contribution_id} is not certified for integration"
            }
        
        # Check if contribution has already been integrated
        if contribution.get("integrated", False):
            return {
                "status": "error",
                "message": f"Contribution {contribution_id} has already been integrated"
            }
        
        # In a real implementation, this would apply the changes to the repository
        # For now, we'll just update the contribution status
        
        # Update contribution status
        contribution["status"] = "integrated"
        contribution["integrated"] = True
        contribution["integration_timestamp"] = datetime.datetime.now().isoformat()
        
        # Log the integration
        self.integration_log["integrations"].append({
            "contribution_id": contribution_id,
            "timestamp": datetime.datetime.now().isoformat(),
            "scale": contribution["scale"],
            "contributor_id": contribution["contributor_id"],
            "value_score": contribution.get("value_metrics", {}).get("overall_score", 1.0)
        })
        
        # Save changes
        self._save_contributions()
        self._save_integration_log()
        
        return {
            "status": "success",
            "message": f"Contribution {contribution_id} integrated successfully",
            "scale": contribution["scale"],
            "value_score": contribution.get("value_metrics", {}).get("overall_score", 1.0)
        }
    
    def update_readme(self):
        """Update the README.md file with information about recent integrations."""
        readme_path = os.path.join(REPO_PATH, "README.md")
        
        # Read current README content
        if os.path.exists(readme_path):
            with open(readme_path, 'r') as f:
                content = f.read()
        else:
            content = "# Synchronism\n\nUnified model bridging quantum mechanics and cosmic evolution through intent dynamics\n\n"
        
        # Get recent integrations
        recent_integrations = sorted(
            self.integration_log["integrations"],
            key=lambda x: x["timestamp"],
            reverse=True
        )[:10]  # Get 10 most recent
        
        # Create integration section
        integration_section = "\n## Recent Contributions\n\n"
        if recent_integrations:
            for integration in recent_integrations:
                date = integration["timestamp"].split("T")[0]
                scale = integration["scale"].capitalize()
                integration_section += f"- [{date}] {scale} scale contribution integrated\n"
        else:
            integration_section += "No recent contributions.\n"
        
        # Add AI-generated line
        integration_section += "\nThis line was added by the AI agent.\n"
        
        # Check if there's already an integration section
        if "## Recent Contributions" in content:
            # Replace existing section
            content = re.sub(
                r"\n## Recent Contributions\n\n(.*?)\n\n",
                f"\n{integration_section}\n\n",
                content,
                flags=re.DOTALL
            )
        else:
            # Add new section
            content += f"\n{integration_section}\n"
        
        # Write updated content
        with open(readme_path, 'w') as f:
            f.write(content)
        
        return {
            "status": "success",
            "message": "README.md updated with recent integrations"
        }
    
    def commit_changes(self):
        """Commit changes to the repository."""
        try:
            # Set Git configuration
            subprocess.run(["git", "config", "--local", "user.email", "actions@github.com"], check=True)
            subprocess.run(["git", "config", "--local", "user.name", "GitHub Actions"], check=True)
            
            # Add changes
            subprocess.run(["git", "add", "."], check=True)
            
            # Create commit message
            commit_message = f"Automated AI Update: Synchronism revision [{datetime.datetime.now().strftime('%Y-%m-%d')}]"
            
            # Commit changes
            result = subprocess.run(
                ["git", "commit", "-m", commit_message],
                capture_output=True,
                text=True
            )
            
            # Check if there were changes to commit
            if "nothing to commit" in result.stdout or "nothing to commit" in result.stderr:
                return {
                    "status": "info",
                    "message": "No changes to commit"
                }
            
            # Push changes
            subprocess.run(["git", "push"], check=True)
            
            return {
                "status": "success",
                "message": "Changes committed and pushed successfully"
            }
            
        except subprocess.CalledProcessError as e:
            return {
                "status": "error",
                "message": f"Git operation failed: {str(e)}",
                "details": e.stderr if hasattr(e, 'stderr') else None
            }


def main():
    """Main function to run when script is executed directly."""
    print("Synchronism Integration System")
    print("=============================")
    
    # Create necessary directories
    os.makedirs(CONFIG_PATH, exist_ok=True)
    
    integration_system = IntegrationSystem()
    
    # Get certified contributions
    certified = integration_system.get_certified_contributions()
    print(f"Found {len(certified)} certified contributions ready for integration")
    
    # Integrate all certified contributions
    for contribution in certified:
        result = integration_system.integrate_contribution(contribution["id"])
        print(f"Integration result: {result['status']} - {result['message']}")
    
    # Update README
    readme_result = integration_system.update_readme()
    print(f"README update: {readme_result['status']} - {readme_result['message']}")
    
    # Commit changes
    commit_result = integration_system.commit_changes()
    print(f"Commit result: {commit_result['status']} - {commit_result['message']}")


if __name__ == "__main__":
    main()
