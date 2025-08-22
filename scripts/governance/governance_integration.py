#!/usr/bin/env python3
"""
Integration module to connect section-specific governance rules 
with the existing governance cycle.
"""

import json
from pathlib import Path
from typing import Dict, Optional, List
from datetime import datetime

from section_governance import SectionGovernance

class GovernanceIntegration:
    """Integrates section-specific rules into the governance workflow."""
    
    def __init__(self, project_root: Path):
        self.project_root = project_root
        self.section_gov = SectionGovernance(project_root)
        
        # Patch existing governance cycle to use section rules
        self._patch_governance_cycle()
    
    def _patch_governance_cycle(self):
        """Prepare integration hooks for governance cycle.
        Note: Actual patching would be done at runtime when cycle is instantiated."""
        # This is a placeholder for now
        # In production, we'd integrate at the cycle creation point
        pass
    
    def evaluate_proposal_with_section_rules(self, proposal: Dict) -> Dict:
        """Enhanced evaluation using section-specific rules."""
        # Get section from proposal
        section = proposal.get('section', '')
        if not section:
            # Try to infer from content
            section = self.analyze_proposal_for_section(proposal)
            proposal['section'] = section
        
        # Use section governance to evaluate
        evaluation = self.section_gov.evaluate_proposal(proposal)
        
        # Update proposal status based on evaluation
        if evaluation['approved']:
            proposal['status'] = 'accepted'
        else:
            proposal['status'] = 'needs_revision'
        
        # Save evaluation metadata
        proposal['evaluation'] = evaluation
        
        return evaluation
    
    def validate_proposal_submission(self, proposal: Dict, 
                                    submitter_tokens: int = 1000) -> tuple:
        """Validate submission against section rules."""
        return self.section_gov.validate_proposal_submission(proposal, submitter_tokens)
    
    def create_governance_rules_document(self) -> str:
        """
        Create a comprehensive governance rules document for participants.
        This gets included in the API context for all governance calls.
        """
        rules_lines = [
            "# Governance Rules (R6 Framework)",
            "",
            "## Overview",
            "The governance system uses section-specific rules that encode different",
            "levels of 'inertia' or filtering for changes. Core principles have high",
            "inertia (resistant to change), while implementation details are fluid.",
            "",
            "## Section-Specific Thresholds",
            ""
        ]
        
        # Add section rules
        for section_name, section_rules in self.section_gov.rules['sections'].items():
            threshold = section_rules.get('change_threshold', 0.5)
            review_days = section_rules.get('review_period_days', 3)
            token_cost = section_rules.get('token_cost', 50)
            description = section_rules.get('description', '')
            
            rules_lines.extend([
                f"### {section_name}",
                f"- **Approval Threshold**: {threshold * 100}%",
                f"- **Review Period**: {review_days} days",
                f"- **Token Cost**: {token_cost}",
                f"- **Description**: {description}",
                ""
            ])
        
        rules_lines.extend([
            "## Proposal Types",
            "",
            "Different types of changes have different requirements:",
            ""
        ])
        
        # Add proposal type modifiers
        for prop_type, type_rules in self.section_gov.rules['proposal_types'].items():
            modifier = type_rules.get('threshold_modifier', 0)
            token_mod = type_rules.get('token_modifier', 1.0)
            
            rules_lines.extend([
                f"### {prop_type.replace('_', ' ').title()}",
                f"- Threshold adjustment: {'+' if modifier >= 0 else ''}{modifier * 100}%",
                f"- Token cost multiplier: {token_mod}x",
            ])
            
            if type_rules.get('fast_track'):
                rules_lines.append("- ⚡ Fast-track available")
            if type_rules.get('requires_quorum'):
                rules_lines.append("- 🏛️ Requires quorum")
            
            rules_lines.append("")
        
        rules_lines.extend([
            "## Resonance Model",
            "",
            "Sections have different 'resonance frequencies' that determine",
            "how quickly they can change:",
            "",
            "- **Maximum** (f=0.05): Foundational principles, very slow change",
            "- **High** (f=0.2): Core concepts, careful evolution",
            "- **Medium** (f=0.5): Standard content, balanced change",
            "- **Low** (f=1.0): Implementation details, rapid iteration",
            "",
            "## Review Weights",
            "",
            "Reviews are weighted by role:",
            "- Arbiter: 40%",
            "- Reviewer: 30%",
            "- Community: 30%",
            "",
            "## Token Economy",
            "",
            "- Base allocation: 1000 tokens per participant",
            "- Successful proposals: Return stake + 20% bonus",
            "- Failed proposals: Lose 50% of stake",
            "- Reviewing: Earn 10-25 tokens based on quality",
            "",
            "## Special Rules",
            "",
            "Some sections have special requirements:",
            "- **Hermetic Principles**: Requires philosophical justification",
            "- **Mathematical Appendix**: Requires mathematical validation",
            "- **Fundamental Changes**: Requires quorum of active participants",
            ""
        ])
        
        return "\n".join(rules_lines)
    
    def update_governance_rules_file(self):
        """Write the governance rules document to disk."""
        rules_file = self.project_root / "scripts" / "governance" / "governance_rules.md"
        rules_content = self.create_governance_rules_document()
        
        with open(rules_file, 'w') as f:
            f.write(rules_content)
        
        print(f"Updated governance rules at: {rules_file}")
        return rules_file
    
    def analyze_proposal_for_section(self, proposal: Dict) -> str:
        """
        Analyze a proposal to determine which section it affects.
        Used when section is not explicitly specified.
        """
        content = proposal.get('content', '').lower()
        title = proposal.get('title', '').lower()
        
        # Keywords for each section
        section_keywords = {
            '03-hermetic-principles': ['hermetic', 'as above', 'correspondence', 'vibration'],
            '04-fundamental-concepts': ['universe grid', 'markov', 'entity', 'interaction'],
            '05-quantum-macro': ['quantum', 'superposition', 'entanglement', 'decoherence'],
            '06-implications': ['implications', 'ethical', 'philosophical'],
            '08-glossary': ['definition', 'glossary', 'terminology'],
            '09-appendix-mathematical': ['mathematical', 'equation', 'proof', 'theorem']
        }
        
        # Score each section
        scores = {}
        for section, keywords in section_keywords.items():
            score = sum(1 for keyword in keywords 
                       if keyword in content or keyword in title)
            if score > 0:
                scores[section] = score
        
        # Return highest scoring section
        if scores:
            return max(scores, key=scores.get)
        
        # Default to fundamental concepts
        return '04-fundamental-concepts'
    
    def generate_status_report(self) -> Dict:
        """Generate a comprehensive governance status report."""
        report = {
            'timestamp': datetime.now().isoformat(),
            'section_configuration': {},
            'active_proposals': [],
            'recent_decisions': [],
            'token_economy': {
                'total_in_circulation': 0,
                'total_staked': 0,
                'recent_transactions': []
            },
            'participation_metrics': {
                'active_participants': 0,
                'proposals_this_week': 0,
                'reviews_this_week': 0
            }
        }
        
        # Add section configuration
        for section in self.section_gov.rules['sections']:
            rules = self.section_gov.get_section_rules(section)
            resonance = self.section_gov.get_resonance_characteristics(section)
            
            report['section_configuration'][section] = {
                'threshold': rules['change_threshold'],
                'resonance_frequency': resonance['frequency'],
                'token_cost': rules['token_cost']
            }
        
        # Would add actual proposal/token data from database
        # This is a template for the full implementation
        
        return report


if __name__ == "__main__":
    # Test the integration
    project_root = Path(__file__).parent.parent.parent
    integration = GovernanceIntegration(project_root)
    
    # Update governance rules file
    rules_file = integration.update_governance_rules_file()
    print(f"\nGovernance rules updated at: {rules_file}")
    
    # Generate status report
    report = integration.generate_status_report()
    print("\nGovernance Status Report:")
    print(json.dumps(report, indent=2))
    
    # Test section analysis
    test_proposal = {
        'title': 'Enhance Markov Blanket Description',
        'content': 'The markov blanket concept in entity interactions needs clarification...'
    }
    
    detected_section = integration.analyze_proposal_for_section(test_proposal)
    print(f"\nProposal would affect section: {detected_section}")
    
    # Show rules for that section
    rules = integration.section_gov.get_section_rules(detected_section)
    print(f"Rules for {detected_section}:")
    print(f"- Threshold: {rules['change_threshold']}")
    print(f"- Review period: {rules['review_period_days']} days")
    print(f"- Token cost: {rules['token_cost']}")