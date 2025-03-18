#!/usr/bin/env python3
"""
Token System Script for Synchronism Governance System

This script manages the ATP/ADP-like token system for the Synchronism repository.
It handles token distribution, exchange, and tracking.
"""

import os
import sys
import json
import datetime
import random
from pathlib import Path

# Constants
# this script is run from Synchronism/scripts/governance so the paths are accordingly
REPO_PATH = os.getenv("REPO_PATH", os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

CONFIG_PATH = os.path.join(REPO_PATH, "scripts", "governance", "config")
TOKENS_PATH = os.path.join(CONFIG_PATH, "tokens.json")
CONTRIBUTIONS_PATH = os.path.join(CONFIG_PATH, "contributions.json")
TENSORS_PATH = os.path.join(CONFIG_PATH, "tensors.json")

# Fractal Scale Tags
FRACTAL_SCALES = {
    "quantum": "Quantum scale phenomena and models",
    "molecular": "Molecular scale interactions and structures",
    "biospheric": "Biological systems and ecosystems",
    "planetary": "Planetary scale processes and systems",
    "galactic": "Cosmic scale phenomena and models"
}

class TokenSystem:
    """Manages the ATP/ADP token system for the Synchronism repository."""
    
    def __init__(self):
        """Initialize the token system."""
        self.tokens = self._load_tokens()
        self.contributions = self._load_contributions()
        self.tensors = self._load_tensors()
        
    def _load_tokens(self):
        """Load token data from file or create if not exists."""
        os.makedirs(CONFIG_PATH, exist_ok=True)
        
        if os.path.exists(TOKENS_PATH):
            with open(TOKENS_PATH, 'r') as f:
                return json.load(f)
        else:
            # Initialize with empty token data
            tokens = {
                "contributors": {},
                "total_supply": {scale: 100 for scale in FRACTAL_SCALES},  # Initial token supply
                "last_distribution": None
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
        """Load contribution data from file."""
        if os.path.exists(CONTRIBUTIONS_PATH):
            with open(CONTRIBUTIONS_PATH, 'r') as f:
                return json.load(f)
        else:
            return {"contributions": []}
    
    def _load_tensors(self):
        """Load tensor data from file."""
        if os.path.exists(TENSORS_PATH):
            with open(TENSORS_PATH, 'r') as f:
                return json.load(f)
        else:
            return {"contributors": {}}
    
    def get_contributor_tokens(self, contributor_id):
        """Get token balance for a contributor."""
        if contributor_id not in self.tokens["contributors"]:
            # Initialize new contributor
            self.tokens["contributors"][contributor_id] = {
                "charged": {scale: 1 for scale in FRACTAL_SCALES},  # Start with 1 token per scale
                "discharged": {scale: 0 for scale in FRACTAL_SCALES},
                "last_activity": datetime.datetime.now().isoformat()
            }
            self._save_tokens()
            
        return self.tokens["contributors"][contributor_id]
    
    def distribute_tokens(self):
        """
        Distribute new tokens to contributors based on their T3/V3 tensor values.
        
        This function is called periodically to ensure token circulation.
        """
        print("Distributing tokens to contributors...")
        
        # Record distribution time
        distribution_time = datetime.datetime.now().isoformat()
        self.tokens["last_distribution"] = distribution_time
        
        # Get all contributors from tensors
        for contributor_id, tensor_data in self.tensors.get("contributors", {}).items():
            # Ensure contributor exists in token system
            tokens = self.get_contributor_tokens(contributor_id)
            
            # Calculate distribution amount based on T3 tensor
            t3 = tensor_data.get("T3", {})
            
            for scale in FRACTAL_SCALES:
                # Calculate scale-specific distribution
                talent = t3.get("talent", {}).get(scale, 0.5)
                training = t3.get("training", {}).get(scale, 0.5)
                temperament = t3.get("temperament", 0.5)
                
                # Weight the components
                trust_score = (0.4 * talent + 0.4 * training + 0.2 * temperament)
                
                # Base distribution amount
                base_amount = 0.5
                
                # Scale by trust score (0.5 to 1.5)
                distribution_amount = base_amount + trust_score
                
                # Add tokens to contributor's balance
                tokens["charged"][scale] = tokens["charged"].get(scale, 0) + distribution_amount
                
                # Update total supply
                self.tokens["total_supply"][scale] = self.tokens["total_supply"].get(scale, 100) + distribution_amount
            
            # Update last activity
            tokens["last_activity"] = distribution_time
        
        # Save changes
        self._save_tokens()
        
        return {
            "status": "success",
            "message": "Tokens distributed successfully",
            "distribution_time": distribution_time,
            "contributor_count": len(self.tensors.get("contributors", {}))
        }
    
    def process_token_exchanges(self):
        """
        Process token exchanges for certified contributions.
        
        This converts discharged tokens back to charged tokens based on value certification.
        """
        print("Processing token exchanges for certified contributions...")
        
        # Find certified contributions that haven't been processed for token exchange
        for contribution in self.contributions["contributions"]:
            if (contribution["status"] == "certified" and 
                contribution.get("value_certified", False) and 
                not contribution.get("token_exchange_processed", False)):
                
                contributor_id = contribution["contributor_id"]
                scale = contribution["scale"]
                value_score = contribution.get("value_metrics", {}).get("overall_score", 1.0)
                
                # Get contributor tokens
                tokens = self.get_contributor_tokens(contributor_id)
                
                # Reduce discharged tokens
                if tokens["discharged"].get(scale, 0) > 0:
                    tokens["discharged"][scale] -= 1
                
                # Add charged tokens based on value score
                # Value score is typically 0.0 to 2.0, with 1.0 being neutral
                tokens["charged"][scale] = tokens["charged"].get(scale, 0) + value_score
                
                # Mark contribution as processed
                contribution["token_exchange_processed"] = True
                contribution["token_exchange_timestamp"] = datetime.datetime.now().isoformat()
                contribution["token_exchange_value"] = value_score
                
                print(f"Processed token exchange for contribution {contribution['id']}")
                print(f"  Scale: {scale}")
                print(f"  Value score: {value_score}")
                print(f"  New charged token balance: {tokens['charged'][scale]}")
        
        # Save changes
        self._save_tokens()
        
        # Save contribution changes
        with open(CONTRIBUTIONS_PATH, 'w') as f:
            json.dump(self.contributions, f, indent=2)
        
        return {
            "status": "success",
            "message": "Token exchanges processed successfully"
        }
    
    def generate_token_report(self):
        """Generate a report of token distribution and usage."""
        report = {
            "timestamp": datetime.datetime.now().isoformat(),
            "total_supply": self.tokens["total_supply"],
            "last_distribution": self.tokens["last_distribution"],
            "contributor_count": len(self.tokens["contributors"]),
            "total_charged_tokens": {scale: 0 for scale in FRACTAL_SCALES},
            "total_discharged_tokens": {scale: 0 for scale in FRACTAL_SCALES},
            "scale_activity": {scale: 0 for scale in FRACTAL_SCALES}
        }
        
        # Calculate totals
        for contributor_id, token_data in self.tokens["contributors"].items():
            for scale in FRACTAL_SCALES:
                report["total_charged_tokens"][scale] += token_data["charged"].get(scale, 0)
                report["total_discharged_tokens"][scale] += token_data["discharged"].get(scale, 0)
        
        # Calculate scale activity from contributions
        for contribution in self.contributions["contributions"]:
            scale = contribution["scale"]
            report["scale_activity"][scale] += 1
        
        # Save report
        report_path = os.path.join(CONFIG_PATH, f"token_report_{datetime.datetime.now().strftime('%Y%m%d')}.json")
        with open(report_path, 'w') as f:
            json.dump(report, f, indent=2)
        
        return {
            "status": "success",
            "message": "Token report generated successfully",
            "report_path": report_path,
            "report": report
        }


def main():
    """Main function to run when script is executed directly."""
    print("Synchronism Token System")
    print("=======================")
    
    # Create necessary directories
    os.makedirs(CONFIG_PATH, exist_ok=True)
    
    token_system = TokenSystem()
    
    # Distribute tokens
    distribution_result = token_system.distribute_tokens()
    print(f"Distribution result: {distribution_result['status']} - {distribution_result['message']}")
    
    # Process token exchanges
    exchange_result = token_system.process_token_exchanges()
    print(f"Exchange result: {exchange_result['status']} - {exchange_result['message']}")
    
    # Generate token report
    report_result = token_system.generate_token_report()
    print(f"Report result: {report_result['status']} - {report_result['message']}")
    print(f"Report saved to: {report_result['report_path']}")


if __name__ == "__main__":
    main()
