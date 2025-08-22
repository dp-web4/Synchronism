#!/usr/bin/env python3
"""Manual review interface for proposals"""

import json
from pathlib import Path
from datetime import datetime

def review_proposal():
    """Interactive review of proposal 003"""
    
    # Load proposals
    config_path = Path("config/whitepaper_proposals.json")
    with open(config_path) as f:
        data = json.load(f)
    
    # Find proposal 003
    proposal = None
    for p in data['proposals']:
        if p['id'] == '003':
            proposal = p
            break
    
    if not proposal:
        print("Proposal 003 not found!")
        return
    
    print("\n" + "="*70)
    print("PROPOSAL FOR REVIEW")
    print("="*70)
    print(f"ID: {proposal['id']}")
    print(f"Title: {proposal['title']}")
    print(f"Author: {proposal['author']}")
    print(f"Type: {proposal['type']}")
    print(f"Section: {proposal['section']}")
    print(f"\nSummary: {proposal['content_summary']}")
    
    # Show full proposal content
    proposal_file = Path("../../whitepaper/sections/04-fundamental-concepts/meta/proposals/003-reorganize-compression-trust.md")
    if proposal_file.exists():
        print("\n" + "="*70)
        print("FULL PROPOSAL CONTENT")
        print("="*70)
        with open(proposal_file) as f:
            print(f.read())
    
    print("\n" + "="*70)
    print("REVIEW OPTIONS")
    print("="*70)
    print("1. ACCEPT - Implement the reorganization as proposed")
    print("2. ACCEPT_WITH_REVISIONS - Agree but suggest modifications")
    print("3. REQUEST_CHANGES - Need significant changes before acceptance")
    print("4. REJECT - Disagree with the proposal")
    print("5. DISCUSS - Continue discussion before decision")
    
    choice = input("\nYour review decision (1-5): ").strip()
    
    recommendation_map = {
        '1': 'accept',
        '2': 'accept_with_revisions', 
        '3': 'request_changes',
        '4': 'reject',
        '5': 'discuss'
    }
    
    if choice not in recommendation_map:
        print("Invalid choice")
        return
    
    recommendation = recommendation_map[choice]
    
    print("\nPlease provide your review comments:")
    print("(Enter 'END' on a new line when finished)")
    
    comments = []
    while True:
        line = input()
        if line == 'END':
            break
        comments.append(line)
    
    review_text = '\n'.join(comments)
    
    # Create review
    review = {
        'reviewer': 'Dennis',
        'date': datetime.now().isoformat(),
        'recommendation': recommendation,
        'comments': review_text,
        'trust_score': 1.0  # Human reviewer has full trust
    }
    
    # Add review to proposal
    proposal['reviews'].append(review)
    proposal['status'] = 'under_review'
    proposal['last_updated'] = datetime.now().isoformat()
    
    # Save updated proposals
    data['last_updated'] = datetime.now().isoformat()
    with open(config_path, 'w') as f:
        json.dump(data, f, indent=2)
    
    print(f"\n✅ Review recorded as: {recommendation}")
    print(f"Comments: {review_text}")
    
    if recommendation in ['accept', 'accept_with_revisions']:
        print("\nAs the arbiter, you can now implement the changes.")
        print("The proposal suggests:")
        print("1. Moving compression-trust content to section 13")
        print("2. Updating the fundamental concepts introduction")
        print("3. Ensuring Universe Grid is presented first")

if __name__ == "__main__":
    review_proposal()