#!/usr/bin/env python3
"""
Contribution Processing Script for Synchronism Governance System

This script handles the processing of new contributions to the Synchronism repository.
It implements the ATP/ADP token mechanism for tracking contributions and their value.
"""

import os
import sys
import json
import hashlib
import datetime
import re
from pathlib import Path

# Constants
# this script is run from Synchronism/scripts/governance so the paths are accordingly
REPO_PATH = os.getenv("REPO_PATH", os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

CONFIG_PATH = os.path.join(REPO_PATH, "scripts", "governance", "config")
TOKENS_PATH = os.path.join(CONFIG_PATH, "tokens.json")
CONTRIBUTIONS_PATH = os.path.join(CONFIG_PATH, "contributions.json")

# Fractal Scale Tags - representing different levels of reality in Synchronism
FRACTAL_SCALES = {
    "quantum": "Quantum scale phenomena and models",
    "molecular": "Molecular scale interactions and structures",
    "biospheric": "Biological systems and ecosystems",
    "planetary": "Planetary scale processes and systems",
    "galactic": "Cosmic scale phenomena and models"
}

class TokenSystem:
    """Manages the ATP/ADP token system for contributions."""
    
    def __init__(self):
        """Initialize the token system."""
        self.tokens = self._load_tokens()
        self.contributions = self._load_contributions()
        
    def _load_tokens(self):
        """Load token data from file or create if not exists."""
        os.makedirs(CONFIG_PATH, exist_ok=True)
        
        if os.path.exists(TOKENS_PATH):
            with open(TOKENS_PATH, 'r') as f:
                return json.load(f)
        else:
            # Initialize with empty token data
            tokens = {
                "contributors": {}
            }
            self._save_tokens(tokens)
            return tokens
    
    def _save_tokens(self, tokens=None):
        """Save token data to file."""
        if tokens is None:
            tokens = self.tokens
            
        with open(TOKENS_PATH, 'w') as f:
            json.dump(tokens, f, indent=2)
    
    def _load_contributions(self):
        """Load contribution data from file or create if not exists."""
        if os.path.exists(CONTRIBUTIONS_PATH):
            with open(CONTRIBUTIONS_PATH, 'r') as f:
                return json.load(f)
        else:
            # Initialize with empty contribution data
            contributions = {
                "contributions": []
            }
            self._save_contributions(contributions)
            return contributions
    
    def _save_contributions(self, contributions=None):
        """Save contribution data to file."""
        if contributions is None:
            contributions = self.contributions
            
        with open(CONTRIBUTIONS_PATH, 'w') as f:
            json.dump(contributions, f, indent=2)
    
    def get_contributor_tokens(self, contributor_id):
        """Get token balance for a contributor."""
        if contributor_id not in self.tokens["contributors"]:
            # Initialize new contributor
            self.tokens["contributors"][contributor_id] = {
                "charged": {scale: 1 for scale in FRACTAL_SCALES},  # Start with 1 token per scale
                "discharged": {scale: 0 for scale in FRACTAL_SCALES}
            }
            self._save_tokens()
            
        return self.tokens["contributors"][contributor_id]
    
    def discharge_token(self, contributor_id, scale, contribution_id):
        """
        Convert a charged token to a discharged token for a contribution.
        
        Args:
            contributor_id: Unique identifier for the contributor
            scale: The fractal scale of the contribution
            contribution_id: Unique identifier for the contribution
            
        Returns:
            bool: True if successful, False if not enough charged tokens
        """
        tokens = self.get_contributor_tokens(contributor_id)
        
        if scale not in tokens["charged"] or tokens["charged"][scale] < 1:
            return False
        
        # Discharge the token
        tokens["charged"][scale] -= 1
        tokens["discharged"][scale] = tokens["discharged"].get(scale, 0) + 1
        
        # Record the contribution
        self.contributions["contributions"].append({
            "id": contribution_id,
            "contributor_id": contributor_id,
            "scale": scale,
            "timestamp": datetime.datetime.now().isoformat(),
            "status": "pending",
            "validations": [],
            "value_certified": False
        })
        
        # Save changes
        self._save_tokens()
        self._save_contributions()
        
        return True
    
    def certify_value(self, contribution_id, validators, value_score):
        """
        Certify the value of a contribution and recharge tokens.
        
        Args:
            contribution_id: Unique identifier for the contribution
            validators: List of validator IDs who certified the value
            value_score: Value score from 0.0 to 2.0 (1.0 is neutral)
            
        Returns:
            bool: True if successful, False if contribution not found
        """
        # Find the contribution
        contribution = None
        for c in self.contributions["contributions"]:
            if c["id"] == contribution_id:
                contribution = c
                break
                
        if contribution is None:
            return False
        
        # Update contribution status
        contribution["status"] = "certified"
        contribution["value_certified"] = True
        contribution["value_score"] = value_score
        contribution["validators"] = validators
        
        # Recharge tokens based on value score
        contributor_id = contribution["contributor_id"]
        scale = contribution["scale"]
        tokens = self.get_contributor_tokens(contributor_id)
        
        # Reduce discharged token
        tokens["discharged"][scale] = max(0, tokens["discharged"][scale] - 1)
        
        # Add charged tokens based on value score
        tokens["charged"][scale] = tokens["charged"].get(scale, 0) + value_score
        
        # Save changes
        self._save_tokens()
        self._save_contributions()
        
        return True


class ContributionProcessor:
    """Processes new contributions to the Synchronism repository."""
    
    def __init__(self):
        """Initialize the contribution processor."""
        self.token_system = TokenSystem()
        
    def generate_contribution_id(self, content, contributor_id):
        """Generate a unique ID for a contribution based on its content and contributor."""
        hash_input = f"{content}{contributor_id}{datetime.datetime.now().isoformat()}"
        return hashlib.sha256(hash_input.encode()).hexdigest()[:16]
    
    def detect_scale(self, content, file_path):
        """
        Detect the fractal scale of a contribution based on content and file path.
        
        Returns the most relevant scale or "quantum" as default.
        """
        # Check file path first
        for scale in FRACTAL_SCALES:
            if scale in file_path.lower():
                return scale
        
        # Check content for scale-related keywords
        scale_keywords = {
            "quantum": ["quantum", "planck", "wave function", "superposition", "entanglement"],
            "molecular": ["molecular", "atom", "chemical", "protein", "dna"],
            "biospheric": ["biosphere", "ecosystem", "organism", "species", "evolution"],
            "planetary": ["planet", "earth", "climate", "geosphere", "atmosphere"],
            "galactic": ["galaxy", "cosmic", "universe", "star", "astronomical"]
        }
        
        scale_counts = {scale: 0 for scale in FRACTAL_SCALES}
        
        for scale, keywords in scale_keywords.items():
            for keyword in keywords:
                scale_counts[scale] += len(re.findall(r'\b' + keyword + r'\b', content, re.IGNORECASE))
        
        # Return the scale with the most keyword matches, or quantum as default
        max_scale = max(scale_counts.items(), key=lambda x: x[1])
        return max_scale[0] if max_scale[1] > 0 else "quantum"
    
    def process_contribution(self, content, file_path, contributor_id):
        """
        Process a new contribution.
        
        Args:
            content: The content of the contribution
            file_path: Path to the file being modified
            contributor_id: Unique identifier for the contributor
            
        Returns:
            dict: Information about the processed contribution
        """
        # Generate a unique ID for this contribution
        contribution_id = self.generate_contribution_id(content, contributor_id)
        
        # Detect the scale of the contribution
        scale = self.detect_scale(content, file_path)
        
        # Discharge a token for this contribution
        success = self.token_system.discharge_token(contributor_id, scale, contribution_id)
        
        if not success:
            return {
                "status": "error",
                "message": f"Not enough charged tokens for scale: {scale}",
                "scale": scale
            }
        
        return {
            "status": "success",
            "contribution_id": contribution_id,
            "scale": scale,
            "message": f"Contribution processed successfully. Token discharged for scale: {scale}"
        }


def main():
    """Main function to run when script is executed directly."""
    print("Synchronism Contribution Processor")
    print("==================================")
    
    # Create necessary directories
    os.makedirs(CONFIG_PATH, exist_ok=True)
    
    processor = ContributionProcessor()
    
    # In a real implementation, this would process contributions from PRs or commits
    # For now, we'll just print a message
    print("Ready to process contributions.")
    print(f"Supported scales: {', '.join(FRACTAL_SCALES.keys())}")
    
    # Example of how this would be used:
    # result = processor.process_contribution(
    #     content="New quantum model for intent transfer",
    #     file_path="Mathematical_Frameworks/Quantum/intent_transfer.md",
    #     contributor_id="github-user123"
    # )
    # print(result)


if __name__ == "__main__":
    main()
