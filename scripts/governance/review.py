#!/usr/bin/env python3
"""
Review System Script for Synchronism Governance System

This script manages the review process for contributions to the Synchronism repository.
It implements a multi-agent review protocol that allows both AI models and human contributors
to review proposed changes.
"""

import os
import sys
import json
import datetime
import random
import uuid
from pathlib import Path

# Constants
# this script is run from Synchronism/scripts/governance so the paths are accordingly
REPO_PATH = os.getenv("REPO_PATH", os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

CONFIG_PATH = os.path.join(REPO_PATH, "scripts", "governance", "config")
REVIEWS_PATH = os.path.join(CONFIG_PATH, "reviews.json")
CONTRIBUTIONS_PATH = os.path.join(CONFIG_PATH, "contributions.json")
TOKENS_PATH = os.path.join(CONFIG_PATH, "tokens.json")

# Fractal Scale Tags
FRACTAL_SCALES = {
    "quantum": "Quantum scale phenomena and models",
    "molecular": "Molecular scale interactions and structures",
    "biospheric": "Biological systems and ecosystems",
    "planetary": "Planetary scale processes and systems",
    "galactic": "Cosmic scale phenomena and models"
}

class ReviewSystem:
    """Manages the review process for contributions to the Synchronism repository."""
    
    def __init__(self):
        """Initialize the review system."""
        self.reviews = self._load_reviews()
        self.contributions = self._load_contributions()
        self.tokens = self._load_tokens()
        
    def _load_reviews(self):
        """Load review data from file or create if not exists."""
        os.makedirs(CONFIG_PATH, exist_ok=True)
        
        if os.path.exists(REVIEWS_PATH):
            with open(REVIEWS_PATH, 'r') as f:
                return json.load(f)
        else:
            # Initialize with empty review data
            reviews = {
                "reviews": [],
                "review_assignments": {},
                "review_metrics": {
                    "total_reviews": 0,
                    "reviews_by_scale": {scale: 0 for scale in FRACTAL_SCALES},
                    "reviews_by_type": {"human": 0, "ai": 0}
                }
            }
            self._save_reviews(reviews)
            return reviews
    
    def _save_reviews(self, reviews=None):
        """Save review data to file."""
        if reviews is None:
            reviews = self.reviews
            
        with open(REVIEWS_PATH, 'w') as f:
            json.dump(reviews, f, indent=2)
    
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
    
    def _load_tokens(self):
        """Load token data from file."""
        if os.path.exists(TOKENS_PATH):
            with open(TOKENS_PATH, 'r') as f:
                return json.load(f)
        else:
            return {"contributors": {}}
    
    def get_pending_reviews(self, reviewer_id=None):
        """
        Get all pending reviews for a specific reviewer or all reviewers.
        
        Args:
            reviewer_id: Optional reviewer ID to filter by
            
        Returns:
            list: Pending review assignments
        """
        assignments = self.reviews.get("review_assignments", {})
        
        if reviewer_id:
            # Get assignments for specific reviewer
            return assignments.get(reviewer_id, [])
        else:
            # Get all assignments
            all_assignments = []
            for reviewer, assignments_list in assignments.items():
                all_assignments.extend(assignments_list)
            return all_assignments
    
    def assign_reviewers(self, contribution_id, num_reviewers=3):
        """
        Assign reviewers to a contribution.
        
        Args:
            contribution_id: Unique identifier for the contribution
            num_reviewers: Number of reviewers to assign
            
        Returns:
            dict: Result of the assignment process
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
        
        # Check if contribution is in the right state for review
        if contribution["status"] != "pending":
            return {
                "status": "error",
                "message": f"Contribution {contribution_id} is not in pending state"
            }
        
        # Get scale of contribution
        scale = contribution["scale"]
        
        # Get potential reviewers (those with tokens in this scale)
        potential_reviewers = []
        for contributor_id, token_data in self.tokens.get("contributors", {}).items():
            # Skip the contributor who made the contribution
            if contributor_id == contribution["contributor_id"]:
                continue
                
            # Check if they have tokens in this scale
            if token_data.get("charged", {}).get(scale, 0) > 0:
                potential_reviewers.append(contributor_id)
        
        # Add AI reviewers
        ai_reviewers = ["ai_reviewer_1", "ai_reviewer_2", "ai_reviewer_3"]
        potential_reviewers.extend(ai_reviewers)
        
        # Shuffle and select reviewers
        random.shuffle(potential_reviewers)
        selected_reviewers = potential_reviewers[:num_reviewers]
        
        # Create review assignments
        for reviewer_id in selected_reviewers:
            # Create assignment
            assignment = {
                "contribution_id": contribution_id,
                "assigned_at": datetime.datetime.now().isoformat(),
                "status": "pending",
                "scale": scale,
                "contributor_id": contribution["contributor_id"]
            }
            
            # Add to reviewer's assignments
            if reviewer_id not in self.reviews["review_assignments"]:
                self.reviews["review_assignments"][reviewer_id] = []
                
            self.reviews["review_assignments"][reviewer_id].append(assignment)
        
        # Update contribution status
        contribution["review_status"] = "assigned"
        contribution["assigned_reviewers"] = selected_reviewers
        contribution["review_assignment_timestamp"] = datetime.datetime.now().isoformat()
        
        # Save changes
        self._save_reviews()
        self._save_contributions()
        
        return {
            "status": "success",
            "message": f"Assigned {len(selected_reviewers)} reviewers to contribution {contribution_id}",
            "reviewers": selected_reviewers
        }
    
    def submit_review(self, reviewer_id, contribution_id, review_data):
        """
        Submit a review for a contribution.
        
        Args:
            reviewer_id: Unique identifier for the reviewer
            contribution_id: Unique identifier for the contribution
            review_data: Dict with review metrics and comments
            
        Returns:
            dict: Result of the review submission
        """
        # Check if reviewer has this assignment
        assignments = self.get_pending_reviews(reviewer_id)
        assignment = None
        
        for a in assignments:
            if a["contribution_id"] == contribution_id:
                assignment = a
                break
                
        if assignment is None:
            return {
                "status": "error",
                "message": f"Reviewer {reviewer_id} does not have an assignment for contribution {contribution_id}"
            }
        
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
        
        # Create review record
        review_id = str(uuid.uuid4())
        review = {
            "id": review_id,
            "reviewer_id": reviewer_id,
            "contribution_id": contribution_id,
            "timestamp": datetime.datetime.now().isoformat(),
            "scale": assignment["scale"],
            "metrics": review_data.get("metrics", {}),
            "comments": review_data.get("comments", ""),
            "recommendation": review_data.get("recommendation", "neutral"),
            "reviewer_type": "ai" if reviewer_id.startswith("ai_") else "human"
        }
        
        # Add review to records
        self.reviews["reviews"].append(review)
        
        # Update review metrics
        self.reviews["review_metrics"]["total_reviews"] += 1
        self.reviews["review_metrics"]["reviews_by_scale"][assignment["scale"]] += 1
        self.reviews["review_metrics"]["reviews_by_type"][review["reviewer_type"]] += 1
        
        # Update assignment status
        assignment["status"] = "completed"
        assignment["review_id"] = review_id
        
        # Check if all reviews are complete
        all_complete = True
        for reviewer in contribution.get("assigned_reviewers", []):
            reviewer_assignments = self.get_pending_reviews(reviewer)
            for a in reviewer_assignments:
                if a["contribution_id"] == contribution_id and a["status"] != "completed":
                    all_complete = False
                    break
        
        # If all reviews are complete, update contribution status
        if all_complete:
            contribution["review_status"] = "completed"
            contribution["review_completion_timestamp"] = datetime.datetime.now().isoformat()
            
            # Calculate review consensus
            self._calculate_review_consensus(contribution)
        
        # Save changes
        self._save_reviews()
        self._save_contributions()
        
        # If reviewer is human, deduct a token
        if review["reviewer_type"] == "human":
            self._deduct_review_token(reviewer_id, assignment["scale"])
        
        return {
            "status": "success",
            "message": f"Review submitted successfully for contribution {contribution_id}",
            "review_id": review_id,
            "all_reviews_complete": all_complete
        }
    
    def _calculate_review_consensus(self, contribution):
        """
        Calculate consensus from reviews for a contribution.
        
        Args:
            contribution: The contribution dict
        """
        # Get all reviews for this contribution
        contribution_reviews = []
        for review in self.reviews["reviews"]:
            if review["contribution_id"] == contribution["id"]:
                contribution_reviews.append(review)
        
        if not contribution_reviews:
            return
        
        # Calculate average metrics
        avg_metrics = {
            "value_score": 0,
            "veracity_score": 0,
            "validity_score": 0,
            "coherence_score": 0,
            "innovation_score": 0
        }
        
        recommendation_counts = {
            "accept": 0,
            "revise": 0,
            "reject": 0,
            "neutral": 0
        }
        
        for review in contribution_reviews:
            metrics = review.get("metrics", {})
            for key in avg_metrics:
                avg_metrics[key] += metrics.get(key, 0.5)
            
            recommendation = review.get("recommendation", "neutral")
            recommendation_counts[recommendation] += 1
        
        # Normalize
        for key in avg_metrics:
            avg_metrics[key] /= len(contribution_reviews)
        
        # Determine consensus recommendation
        max_count = 0
        consensus = "neutral"
        
        for rec, count in recommendation_counts.items():
            if count > max_count:
                max_count = count
                consensus = rec
        
        # Update contribution with consensus
        contribution["review_consensus"] = {
            "metrics": avg_metrics,
            "recommendation": consensus,
            "review_count": len(contribution_reviews),
            "recommendation_counts": recommendation_counts
        }
        
        # Update status based on consensus
        if consensus == "accept":
            contribution["status"] = "accepted"
        elif consensus == "reject":
            contribution["status"] = "rejected"
        elif consensus == "revise":
            contribution["status"] = "needs_revision"
        else:
            # If neutral, use metric thresholds
            overall_score = (avg_metrics["value_score"] + 
                            avg_metrics["veracity_score"] + 
                            avg_metrics["validity_score"]) / 3
            
            if overall_score >= 0.7:
                contribution["status"] = "accepted"
            elif overall_score <= 0.3:
                contribution["status"] = "rejected"
            else:
                contribution["status"] = "needs_revision"
    
    def _deduct_review_token(self, reviewer_id, scale):
        """
        Deduct a token from a reviewer for completing a review.
        
        Args:
            reviewer_id: Unique identifier for the reviewer
            scale: Scale of the contribution reviewed
        """
        # Get reviewer tokens
        if reviewer_id in self.tokens["contributors"]:
            tokens = self.tokens["contributors"][reviewer_id]
            
            # Deduct charged token if available
            if tokens.get("charged", {}).get(scale, 0) > 0:
                tokens["charged"][scale] -= 1
                
                # Add to discharged tokens
                if "discharged" not in tokens:
                    tokens["discharged"] = {}
                if scale not in tokens["discharged"]:
                    tokens["discharged"][scale] = 0
                    
                tokens["discharged"][scale] += 1
                
                # Save changes
                with open(TOKENS_PATH, 'w') as f:
                    json.dump(self.tokens, f, indent=2)
    
    def process_accepted_contributions(self):
        """
        Process contributions that have been accepted through review.
        
        This moves them to the validation stage.
        """
        processed_count = 0
        
        for contribution in self.contributions["contributions"]:
            if contribution["status"] == "accepted" and not contribution.get("processed_for_validation", False):
                # Mark as ready for validation
                contribution["status"] = "pending_validation"
                contribution["validation_ready_timestamp"] = datetime.datetime.now().isoformat()
                contribution["processed_for_validation"] = True
                
                processed_count += 1
        
        # Save changes
        self._save_contributions()
        
        return {
            "status": "success",
            "message": f"Processed {processed_count} accepted contributions for validation"
        }
    
    def generate_review_report(self):
        """Generate a report of review activity."""
        report = {
            "timestamp": datetime.datetime.now().isoformat(),
            "total_reviews": self.reviews["review_metrics"]["total_reviews"],
            "reviews_by_scale": self.reviews["review_metrics"]["reviews_by_scale"],
            "reviews_by_type": self.reviews["review_metrics"]["reviews_by_type"],
            "pending_reviews": len(self.get_pending_reviews()),
            "review_recommendations": {
                "accept": 0,
                "revise": 0,
                "reject": 0,
                "neutral": 0
            }
        }
        
        # Count recommendations
        for review in self.reviews["reviews"]:
            recommendation = review.get("recommendation", "neutral")
            report["review_recommendations"][recommendation] += 1
        
        # Save report
        report_path = os.path.join(CONFIG_PATH, f"review_report_{datetime.datetime.now().strftime('%Y%m%d')}.json")
        with open(report_path, 'w') as f:
            json.dump(report, f, indent=2)
        
        return {
            "status": "success",
            "message": "Review report generated successfully",
            "report_path": report_path,
            "report": report
        }


def main():
    """Main function to run when script is executed directly."""
    print("Synchronism Review System")
    print("========================")
    
    # Create necessary directories
    os.makedirs(CONFIG_PATH, exist_ok=True)
    
    review_system = ReviewSystem()
    
    # Process accepted contributions
    process_result = review_system.process_accepted_contributions()
    print(f"Processing result: {process_result['status']} - {process_result['message']}")
    
    # Generate review report
    report_result = review_system.generate_review_report()
    print(f"Report result: {report_result['status']} - {report_result['message']}")
    print(f"Report saved to: {report_result['report_path']}")
    
    # In a real implementation, this would also assign reviewers to new contributions
    # and process AI reviews automatically


if __name__ == "__main__":
    main()
