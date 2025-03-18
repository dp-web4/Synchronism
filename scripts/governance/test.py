#!/usr/bin/env python3
"""
Test Script for Synchronism Governance System

This script runs tests on the self-governing repository management system
to ensure all components are functioning correctly.
"""

import os
import sys
import json
import datetime
import unittest
from pathlib import Path

# Add the governance scripts directory to the Python path
REPO_PATH = os.getenv("REPO_PATH", os.getcwd())
SCRIPTS_PATH = os.path.join(REPO_PATH, "scripts", "governance")
sys.path.append(SCRIPTS_PATH)

# Import governance modules
try:
    from contribution import ContributionSystem
    from validation import ValidationSystem
    from integration import IntegrationSystem
    from token_system import TokenSystem
    from review import ReviewSystem
    from fractal_branches import FractalBranchSystem
except ImportError as e:
    print(f"Error importing governance modules: {e}")
    print(f"Make sure all required scripts are in {SCRIPTS_PATH}")
    sys.exit(1)

class TestContributionSystem(unittest.TestCase):
    """Test cases for the ContributionSystem."""
    
    def setUp(self):
        """Set up test environment."""
        self.contribution_system = ContributionSystem()
    
    def test_create_contribution(self):
        """Test creating a contribution."""
        test_contribution = {
            "title": "Test Contribution",
            "description": "This is a test contribution for unit testing",
            "content": "Test content for the Synchronism model",
            "scale": "quantum",
            "contributor_id": "test_contributor_unit"
        }
        
        result = self.contribution_system.create_contribution(test_contribution)
        self.assertEqual(result["status"], "success")
        self.assertIn("contribution_id", result)
    
    def test_get_pending_contributions(self):
        """Test getting pending contributions."""
        pending = self.contribution_system.get_pending_contributions()
        self.assertIsInstance(pending, list)


class TestValidationSystem(unittest.TestCase):
    """Test cases for the ValidationSystem."""
    
    def setUp(self):
        """Set up test environment."""
        self.validation_system = ValidationSystem()
        self.contribution_system = ContributionSystem()
        
        # Create a test contribution if needed
        pending = self.contribution_system.get_pending_contributions()
        if not pending:
            test_contribution = {
                "title": "Test Contribution for Validation",
                "description": "This is a test contribution for validation testing",
                "content": "Test content for validation",
                "scale": "molecular",
                "contributor_id": "test_contributor_validation"
            }
            result = self.contribution_system.create_contribution(test_contribution)
            self.contribution_id = result["contribution_id"]
        else:
            self.contribution_id = pending[0]["id"]
    
    def test_validate_contribution(self):
        """Test validating a contribution."""
        validation_data = {
            "value_score": 0.8,
            "veracity_score": 0.7,
            "validity_score": 0.9,
            "comments": "Test validation for unit testing"
        }
        
        result = self.validation_system.validate_contribution(
            contribution_id=self.contribution_id,
            validator_id="test_validator_unit",
            validation_data=validation_data
        )
        
        self.assertEqual(result["status"], "success")


class TestIntegrationSystem(unittest.TestCase):
    """Test cases for the IntegrationSystem."""
    
    def setUp(self):
        """Set up test environment."""
        self.integration_system = IntegrationSystem()
    
    def test_update_readme(self):
        """Test updating the README file."""
        result = self.integration_system.update_readme()
        self.assertEqual(result["status"], "success")


class TestTokenSystem(unittest.TestCase):
    """Test cases for the TokenSystem."""
    
    def setUp(self):
        """Set up test environment."""
        self.token_system = TokenSystem()
    
    def test_get_contributor_tokens(self):
        """Test getting contributor tokens."""
        tokens = self.token_system.get_contributor_tokens("test_contributor_token")
        self.assertIn("charged", tokens)
        self.assertIn("discharged", tokens)
    
    def test_distribute_tokens(self):
        """Test distributing tokens."""
        result = self.token_system.distribute_tokens()
        self.assertEqual(result["status"], "success")


class TestReviewSystem(unittest.TestCase):
    """Test cases for the ReviewSystem."""
    
    def setUp(self):
        """Set up test environment."""
        self.review_system = ReviewSystem()
        self.contribution_system = ContributionSystem()
        
        # Create a test contribution if needed
        pending = self.contribution_system.get_pending_contributions()
        if not pending:
            test_contribution = {
                "title": "Test Contribution for Review",
                "description": "This is a test contribution for review testing",
                "content": "Test content for review",
                "scale": "biospheric",
                "contributor_id": "test_contributor_review"
            }
            result = self.contribution_system.create_contribution(test_contribution)
            self.contribution_id = result["contribution_id"]
        else:
            self.contribution_id = pending[0]["id"]
    
    def test_assign_reviewers(self):
        """Test assigning reviewers to a contribution."""
        result = self.review_system.assign_reviewers(self.contribution_id, num_reviewers=2)
        self.assertEqual(result["status"], "success")


class TestFractalBranchSystem(unittest.TestCase):
    """Test cases for the FractalBranchSystem."""
    
    def setUp(self):
        """Set up test environment."""
        self.branch_system = FractalBranchSystem()
    
    def test_get_branch(self):
        """Test getting branch data."""
        branch = self.branch_system.get_branch("main")
        self.assertIsNotNone(branch)
        self.assertEqual(branch["name"], "main")
    
    def test_generate_branch_report(self):
        """Test generating a branch report."""
        result = self.branch_system.generate_branch_report()
        self.assertEqual(result["status"], "success")


def run_tests():
    """Run all test cases."""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add test cases
    suite.addTests(loader.loadTestsFromTestCase(TestContributionSystem))
    suite.addTests(loader.loadTestsFromTestCase(TestValidationSystem))
    suite.addTests(loader.loadTestsFromTestCase(TestIntegrationSystem))
    suite.addTests(loader.loadTestsFromTestCase(TestTokenSystem))
    suite.addTests(loader.loadTestsFromTestCase(TestReviewSystem))
    suite.addTests(loader.loadTestsFromTestCase(TestFractalBranchSystem))
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return result


if __name__ == "__main__":
    print("Running tests for Synchronism Governance System...")
    result = run_tests()
    
    # Print summary
    print("\nTest Summary:")
    print(f"Ran {result.testsRun} tests")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    
    # Exit with appropriate code
    sys.exit(len(result.failures) + len(result.errors))
