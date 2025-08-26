#!/usr/bin/env python3
"""
Main Script for Synchronism Governance System

This script integrates all components of the self-governing repository management system
and provides a command-line interface for interacting with the system.
"""

import os
import sys
import json
import datetime
import argparse
import subprocess
from pathlib import Path

# Add the governance scripts directory to the Python path
# this script is run from Synchronism/scripts/governance so the paths are accordingly
REPO_PATH = os.getenv("REPO_PATH", os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
SCRIPTS_PATH = os.path.join(REPO_PATH, "scripts", "governance")
sys.path.append(SCRIPTS_PATH)

# Import governance modules
try:
    from contribution import ContributionProcessor
    from validation import ValidationSystem
    from integration import IntegrationSystem
    from token_system import TokenSystem
    from review import ReviewSystem
    from fractal_branches import FractalBranchSystem
    from whitepaper_builder import WhitepaperBuilder
    # Import AI participant modules
    from lct_participants import LCT, LCTRegistry, ParticipantType, AccessMethod
    from participant_api import ClaudeAPI, GPTAPI
    from whitepaper_governance import WhitepaperGovernance
except ImportError as e:
    print(f"Error importing governance modules: {e}")
    print(f"Make sure all required scripts are in {SCRIPTS_PATH}")
    sys.exit(1)

# Constants
CONFIG_PATH = os.path.join(SCRIPTS_PATH, "config")
LOG_PATH = os.path.join(CONFIG_PATH, "governance_log.json")

# Fractal Scale Tags
FRACTAL_SCALES = {
    "quantum": "Quantum scale phenomena and models",
    "molecular": "Molecular scale interactions and structures",
    "biospheric": "Biological systems and ecosystems",
    "planetary": "Planetary scale processes and systems",
    "galactic": "Cosmic scale phenomena and models"
}

class SynchronismGovernanceSystem:
    """Main class for the Synchronism Governance System."""
    
    def __init__(self):
        """Initialize the governance system."""
        # Create necessary directories
        os.makedirs(CONFIG_PATH, exist_ok=True)
        
        # Initialize subsystems in order: token, contribution, validation,
        # integration, review, branch
        self.token_system = TokenSystem()  #initialize first!
        self.contribution_system = ContributionProcessor()
        self.validation_system = ValidationSystem()
        self.integration_system = IntegrationSystem()
        self.review_system = ReviewSystem()
        self.branch_system = FractalBranchSystem()
        
        # Initialize whitepaper governance
        self.whitepaper_governance = WhitepaperGovernance()
        
        # Initialize LCT registry
        self.lct_registry = LCTRegistry()
        
        # Initialize log
        self.log = self._load_log()
    
    def _load_log(self):
        """Load governance log from file or create if not exists."""
        if os.path.exists(LOG_PATH):
            with open(LOG_PATH, 'r') as f:
                return json.load(f)
        else:
            # Initialize with empty log
            log = {
                "events": [],
                "last_run": None,
                "system_status": "initialized",
                "version": "1.0.0"
            }
            self._save_log(log)
            return log
    
    def _save_log(self, log=None):
        """Save governance log to file."""
        if log is None:
            log = self.log
            
        with open(LOG_PATH, 'w') as f:
            json.dump(log, f, indent=2)
    
    def _log_event(self, event_type, details):
        """Log an event to the governance log."""
        event = {
            "timestamp": datetime.datetime.now().isoformat(),
            "type": event_type,
            "details": details
        }
        
        self.log["events"].append(event)
        self.log["last_run"] = event["timestamp"]
        
        self._save_log()
    
    def _get_ai_participants(self):
        """Initialize and return AI participants for proposal generation."""
        participants = []
        
        # Check for Claude API key
        claude_key = os.getenv("ANTHROPIC_API_KEY")
        if claude_key:
            claude_lct = self.lct_registry.get_or_create(
                name="Claude-4.1",
                participant_type=ParticipantType.AI_CLAUDE,
                access_method=AccessMethod.API,
                access_endpoint="https://api.anthropic.com",
                access_credentials={"api_key": claude_key}
            )
            participants.append(ClaudeAPI(claude_lct))
            print("  ✓ Claude participant initialized")
        else:
            print("  ⚠ Claude API key not found - skipping Claude participant")
        
        # Check for GPT API key
        gpt_key = os.getenv("OPENAI_API_KEY")
        if gpt_key:
            gpt_lct = self.lct_registry.get_or_create(
                name="GPT-5",
                participant_type=ParticipantType.AI_GPT,
                access_method=AccessMethod.API,
                access_endpoint="https://api.openai.com",
                access_credentials={"api_key": gpt_key}
            )
            participants.append(GPTAPI(gpt_lct))
            print("  ✓ GPT participant initialized")
        else:
            print("  ⚠ GPT API key not found - skipping GPT participant")
        
        return participants
    
    def run_daily_update(self):
        """
        Run the daily update process for the governance system.
        
        This function orchestrates the entire governance workflow:
        0. AI participants generate proposals
        1. Process new contributions
        2. Assign reviewers to pending contributions
        3. Process reviews and determine consensus
        4. Validate accepted contributions
        5. Integrate validated contributions
        6. Distribute tokens based on contribution value
        7. Update repository with changes
        
        Returns:
            dict: Result of the daily update process
        """
        print("Running daily update for Synchronism Governance System...")
        
        # Start time
        start_time = datetime.datetime.now()
        
        # Log start of daily update
        self._log_event("daily_update_start", {
            "timestamp": start_time.isoformat()
        })
        
        # 0. AI Participants generate proposals
        print("\n0. AI Participants generating proposals...")
        proposals_created = []
        ai_participants = self._get_ai_participants()
        
        if ai_participants:
            # Define sections that can be analyzed for proposals
            whitepaper_sections = [
                "01-hermetic-principles",
                "02-introductory-summary",
                "03-historical-evolution",
                "04-fundamental-concepts",
                "05-mathematical-frameworks",
                "06-universe-as-consciousness-network",
                "07-scale-specific-implementations",
                "08-emergent-properties",
                "09-governance-structure",
                "glossary"
            ]
            
            # Each participant generates up to 1 proposal
            max_proposals_per_participant = 1
            
            for participant in ai_participants:
                try:
                    print(f"\n  {participant.lct.name} analyzing sections...")
                    # Generate proposals (the API will handle the actual analysis)
                    proposals = participant.generate_proposals(
                        sections=whitepaper_sections,
                        max_proposals=max_proposals_per_participant,
                        context={
                            "cycle_id": self.log.get("cycle_count", 0) + 1,
                            "existing_proposals": proposals_created
                        }
                    )
                    
                    if proposals:
                        proposals_created.extend(proposals)
                        print(f"  {participant.lct.name} created {len(proposals)} proposal(s)")
                        for proposal in proposals:
                            print(f"    - {proposal.get('title', 'Untitled')}")
                    else:
                        print(f"  {participant.lct.name} chose not to create any proposals")
                    
                    # Update participant activity
                    participant.update_activity()
                    
                except Exception as e:
                    print(f"  ⚠ Error generating proposals for {participant.lct.name}: {str(e)}")
        else:
            print("  ⚠ No AI participants available for proposal generation")
        
        print(f"\nTotal proposals created: {len(proposals_created)}")
        
        # 1. Process new contributions
        print("\n1. Processing new contributions...")
        pending_contributions = self.validation_system.get_pending_contributions()
        contribution_result = self.contribution_system.process_new_contributions(pending_contributions)
        print(f"Contribution result: {contribution_result['status']} - {contribution_result['message']}")
        
        # 2. Assign reviewers to pending contributions
        print("\n2. Assigning reviewers to pending contributions...")
        pending_contributions = self.validation_system.get_pending_contributions()
        for contribution in pending_contributions:
            if not contribution.get("review_status") == "assigned":
                review_result = self.review_system.assign_reviewers(contribution["id"])
                print(f"Review assignment result: {review_result['status']} - {review_result['message']}")
        
        # 3. Process reviews and determine consensus
        print("\n3. Processing reviews and determining consensus...")
        process_result = self.review_system.process_accepted_contributions()
        print(f"Review processing result: {process_result['status']} - {process_result['message']}")
        
        # 4. Validate accepted contributions
        print("\n4. Validating accepted contributions...")
        pending_validations = self.validation_system.get_pending_contributions()
        print(f"Found {len(pending_validations)} contributions pending validation")
        
        # 5. Integrate validated contributions
        print("\n5. Integrating validated contributions...")
        certified_contributions = self.integration_system.get_certified_contributions()
        for contribution in certified_contributions:
            integration_result = self.integration_system.integrate_contribution(contribution["id"])
            print(f"Integration result: {integration_result['status']} - {integration_result['message']}")
        
        # 6. Distribute tokens based on contribution value
        print("\n6. Distributing tokens based on contribution value...")
        token_result = self.token_system.distribute_tokens()
        print(f"Token distribution result: {token_result['status']} - {token_result['message']}")
        
        exchange_result = self.token_system.process_token_exchanges()
        print(f"Token exchange result: {exchange_result['status']} - {exchange_result['message']}")
        
        # 7. Update repository with changes
        print("\n7. Updating repository with changes...")
        readme_result = self.integration_system.update_readme()
        print(f"README update result: {readme_result['status']} - {readme_result['message']}")
        
        commit_result = self.integration_system.commit_changes()
        print(f"Commit result: {commit_result['status']} - {commit_result['message']}")
        
        # 8. Check for whitepaper rebuilds if fractal files changed
        print("\n8. Checking for whitepaper rebuild...")
        builder = WhitepaperBuilder(REPO_PATH)
        build_result = builder.rebuild_if_changed()
        
        if build_result:
            print(f"Whitepaper rebuilt: {build_result['timestamp']}")
            # Log successful formats
            for fmt, info in build_result.get('formats', {}).items():
                if info['success']:
                    print(f"  ✅ {fmt}: {info['message']}")
        else:
            print("No fractal file changes - whitepaper rebuild skipped")
        
        # End time
        end_time = datetime.datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        # Generate reports
        print("\nGenerating reports...")
        token_report = self.token_system.generate_token_report()
        review_report = self.review_system.generate_review_report()
        branch_report = self.branch_system.generate_branch_report()
        
        # Log completion of daily update
        self._log_event("daily_update_complete", {
            "start_time": start_time.isoformat(),
            "end_time": end_time.isoformat(),
            "duration_seconds": duration,
            "proposals_generated": len(proposals_created),
            "contributions_processed": len(pending_contributions),
            "validations_processed": len(pending_validations),
            "integrations_processed": len(certified_contributions),
            "whitepaper_rebuilt": build_result is not None,
            "ai_participants": len(ai_participants) if ai_participants else 0
        })
        
        return {
            "status": "success",
            "message": "Daily update completed successfully",
            "duration_seconds": duration,
            "reports": {
                "token_report": token_report["report_path"],
                "review_report": review_report["report_path"],
                "branch_report": branch_report["report_path"]
            }
        }
    
    def initialize_system(self):
        """
        Initialize the governance system.
        
        This function sets up the initial state of the system:
        1. Create necessary directories and configuration files
        2. Initialize the fractal branch structure
        3. Set up initial token distribution
        
        Returns:
            dict: Result of the initialization process
        """
        print("Initializing Synchronism Governance System...")
        
        # Log initialization
        self._log_event("system_initialization", {
            "timestamp": datetime.datetime.now().isoformat()
        })
        
        # Initialize fractal branch structure
        print("\nInitializing fractal branch structure...")
        branch_result = self.branch_system.initialize_fractal_structure()
        print(f"Branch initialization result: {branch_result['status']} - {branch_result['message']}")
        
        # Initialize token distribution
        print("\nInitializing token distribution...")
        token_result = self.token_system.distribute_tokens()
        print(f"Token initialization result: {token_result['status']} - {token_result['message']}")
        
        # Update system status
        self.log["system_status"] = "ready"
        self._save_log()
        
        return {
            "status": "success",
            "message": "Synchronism Governance System initialized successfully",
            "branch_result": branch_result,
            "token_result": token_result
        }
    
    def test_system(self):
        """
        Test the governance system functionality.
        
        This function runs tests on each component of the system:
        1. Test contribution submission
        2. Test review assignment and processing
        3. Test validation
        4. Test integration
        5. Test token distribution
        6. Test branch management
        
        Returns:
            dict: Result of the test process
        """
        print("Testing Synchronism Governance System...")
        
        # Log test start
        self._log_event("system_test_start", {
            "timestamp": datetime.datetime.now().isoformat()
        })
        
        test_results = {
            "contribution": self._test_contribution(),
            "review": self._test_review(),
            "validation": self._test_validation(),
            "integration": self._test_integration(),
            "token": self._test_token(),
            "branch": self._test_branch()
        }
        
        # Calculate overall test result
        success_count = sum(1 for result in test_results.values() if result["status"] == "success")
        total_count = len(test_results)
        
        overall_status = "success" if success_count == total_count else "partial" if success_count > 0 else "failure"
        
        # Log test completion
        self._log_event("system_test_complete", {
            "timestamp": datetime.datetime.now().isoformat(),
            "overall_status": overall_status,
            "success_count": success_count,
            "total_count": total_count,
            "results": test_results
        })
        
        return {
            "status": overall_status,
            "message": f"System test completed with {success_count}/{total_count} successful components",
            "results": test_results
        }
    
    def _test_contribution(self):
        """Test the contribution system."""
        print("\nTesting contribution system...")
        try:
            # Test creating a contribution
            test_contribution = {
                "title": "Test Contribution",
                "description": "This is a test contribution for system testing",
                "content": "Test content for the Synchronism model",
                "scale": "quantum",
                "contributor_id": "test_contributor"
            }
            
            result = self.contribution_system.create_contribution(test_contribution)
            
            return {
                "status": "success",
                "message": "Contribution system test passed",
                "details": result
            }
        except Exception as e:
            return {
                "status": "failure",
                "message": f"Contribution system test failed: {str(e)}"
            }
    
    def _test_review(self):
        """Test the review system."""
        print("\nTesting review system...")
        try:
            # Get pending contributions
            pending = self.validation_system.get_pending_contributions()
            
            if not pending:
                return {
                    "status": "skipped",
                    "message": "No pending contributions to test review system"
                }
            
            # Assign reviewers to first pending contribution
            contribution_id = pending[0]["id"]
            result = self.review_system.assign_reviewers(contribution_id, num_reviewers=2)
            
            return {
                "status": "success",
                "message": "Review system test passed",
                "details": result
            }
        except Exception as e:
            return {
                "status": "failure",
                "message": f"Review system test failed: {str(e)}"
            }
    
    def _test_validation(self):
        """Test the validation system."""
        print("\nTesting validation system...")
        try:
            # Create a test validation
            validation_data = {
                "value_score": 0.8,
                "veracity_score": 0.7,
                "validity_score": 0.9,
                "comments": "Test validation for system testing"
            }
            
            # Get pending contributions
            pending = self.validation_system.get_pending_contributions()
            
            if not pending:
                return {
                    "status": "skipped",
                    "message": "No pending contributions to test validation system"
                }
            
            # Validate first pending contribution
            contribution_id = pending[0]["id"]
            result = self.validation_system.validate_contribution(
                contribution_id=contribution_id,
                validator_id="test_validator",
                validation_data=validation_data
            )
            
            return {
                "status": "success",
                "message": "Validation system test passed",
                "details": result
            }
        except Exception as e:
            return {
                "status": "failure",
                "message": f"Validation system test failed: {str(e)}"
            }
    
    def _test_integration(self):
        """Test the integration system."""
        print("\nTesting integration system...")
        try:
            # Test updating README
            result = self.integration_system.update_readme()
            
            return {
                "status": "success",
                "message": "Integration system test passed",
                "details": result
            }
        except Exception as e:
            return {
                "status": "failure",
                "message": f"Integration system test failed: {str(e)}"
            }
    
    def _test_token(self):
        """Test the token system."""
        print("\nTesting token system...")
        try:
            # Test token distribution
            result = self.token_system.distribute_tokens()
            
            return {
                "status": "success",
                "message": "Token system test passed",
                "details": result
            }
        except Exception as e:
            return {
                "status": "failure",
                "message": f"Token system test failed: {str(e)}"
            }
    
    def _test_branch(self):
        """Test the branch system."""
        print("\nTesting branch system...")
        try:
            # Test generating branch report
            result = self.branch_system.generate_branch_report()
            
            return {
                "status": "success",
                "message": "Branch system test passed",
                "details": result
            }
        except Exception as e:
            return {
                "status": "failure",
                "message": f"Branch system test failed: {str(e)}"
            }


def main():
    """Main function to run when script is executed directly."""
    parser = argparse.ArgumentParser(description="Synchronism Governance System")
    parser.add_argument("command", choices=["init", "update", "test"], 
                        help="Command to execute: init (initialize system), update (run daily update), test (test system)")
    
    args = parser.parse_args()
    
    # Create governance system
    governance_system = SynchronismGovernanceSystem()
    
    # Execute command
    if args.command == "init":
        result = governance_system.initialize_system()
    elif args.command == "update":
        result = governance_system.run_daily_update()
    elif args.command == "test":
        result = governance_system.test_system()
    
    print(f"\nCommand '{args.command}' completed with status: {result['status']}")
    print(f"Message: {result['message']}")


if __name__ == "__main__":
    main()
