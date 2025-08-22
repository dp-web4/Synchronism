#!/usr/bin/env python3
"""
Section-specific governance rules implementation.
Implements the R6 workflow with section-specific inertia/filtering.
"""

import json
from pathlib import Path
from typing import Dict, Optional, Tuple, List
from datetime import datetime, timedelta
import logging

class SectionGovernance:
    """Manages section-specific governance rules and thresholds."""
    
    def __init__(self, project_root: Path):
        self.project_root = project_root
        self.rules_file = project_root / "governance" / "section_rules.json"
        self.proposals_dir = project_root / "governance" / "proposals"
        self.rules = self._load_rules()
        
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)
    
    def _load_rules(self) -> Dict:
        """Load governance rules from JSON."""
        if not self.rules_file.exists():
            raise FileNotFoundError(f"Rules file not found: {self.rules_file}")
        
        with open(self.rules_file, 'r') as f:
            return json.load(f)['governance_rules']
    
    def get_section_rules(self, section_path: str) -> Dict:
        """
        Get rules for a specific section.
        
        Args:
            section_path: Path like "04-fundamental-concepts/01-universe-grid"
        
        Returns:
            Dict containing applicable rules with inheritance
        """
        # Start with default rules
        rules = self.rules['default'].copy()
        
        # Parse section path
        parts = section_path.split('/')
        main_section = parts[0] if parts else None
        subsection = parts[1] if len(parts) > 1 else None
        
        # Apply main section rules
        if main_section in self.rules['sections']:
            section_rules = self.rules['sections'][main_section]
            rules.update({k: v for k, v in section_rules.items() 
                         if k != 'subsection_rules'})
            
            # Apply subsection rules if they exist
            if subsection and 'subsection_rules' in section_rules:
                if subsection in section_rules['subsection_rules']:
                    rules.update(section_rules['subsection_rules'][subsection])
        
        return rules
    
    def calculate_approval_threshold(self, 
                                    section_path: str,
                                    proposal_type: str = "addition") -> Tuple[float, Dict]:
        """
        Calculate the approval threshold for a proposal.
        
        Args:
            section_path: Section being modified
            proposal_type: Type of change (typo_fix, clarification, etc.)
        
        Returns:
            Tuple of (threshold, full_rules_dict)
        """
        # Get base section rules
        rules = self.get_section_rules(section_path)
        base_threshold = rules['change_threshold']
        
        # Apply proposal type modifiers
        if proposal_type in self.rules['proposal_types']:
            type_rules = self.rules['proposal_types'][proposal_type]
            threshold = base_threshold + type_rules.get('threshold_modifier', 0)
            threshold = max(0.1, min(1.0, threshold))  # Clamp to [0.1, 1.0]
            
            # Add type-specific rules
            rules['proposal_type'] = proposal_type
            rules['token_cost'] = int(rules['token_cost'] * 
                                     type_rules.get('token_modifier', 1.0))
            rules['fast_track'] = type_rules.get('fast_track', False)
            rules['requires_quorum'] = type_rules.get('requires_quorum', False)
        else:
            threshold = base_threshold
        
        rules['calculated_threshold'] = threshold
        return threshold, rules
    
    def get_resonance_characteristics(self, section_path: str) -> Dict:
        """
        Get the full LRC resonance characteristics for a section.
        
        Returns dict with:
        - L (inductance): Resists change
        - C (capacitance): Stores potential for change
        - R (resistance): Dissipates bad proposals
        - Q_factor: Quality factor (sharpness of resonance)
        - damping: System damping characteristic
        - rejection_rate: How aggressively bad proposals are filtered
        """
        rules = self.get_section_rules(section_path)
        resonance_level = rules.get('resonance', 'medium')
        
        if resonance_level in self.rules['resonance_bands']:
            return self.rules['resonance_bands'][resonance_level]
        
        # Default to medium resonance
        return self.rules['resonance_bands']['medium']
    
    def calculate_proposal_energy_dissipation(self, proposal: Dict, 
                                            rejection_type: str = "review_rejection") -> Dict:
        """
        Calculate energy dissipation for a rejected proposal.
        
        Energy dissipation prevents the system from accumulating bad proposals
        and ensures stability through the resistance mechanism.
        
        Args:
            proposal: The proposal being rejected
            rejection_type: Type of rejection mechanism
        
        Returns:
            Dict with energy_loss and token_return percentages
        """
        section = proposal.get('section', '')
        resonance = self.get_resonance_characteristics(section)
        
        # Get base rejection mechanism
        if rejection_type in self.rules.get('rejection_mechanisms', {}):
            mechanism = self.rules['rejection_mechanisms'][rejection_type]
        else:
            mechanism = self.rules['rejection_mechanisms']['review_rejection']
        
        # Adjust based on section's R value
        r_value = resonance.get('R', 2)
        rejection_rate = resonance.get('rejection_rate', 0.3)
        
        # Higher R = more dissipation for that section
        adjusted_energy_loss = mechanism['energy_loss'] * (1 + rejection_rate)
        adjusted_energy_loss = min(1.0, adjusted_energy_loss)  # Cap at 100%
        
        return {
            'energy_loss': adjusted_energy_loss,
            'token_return': 1.0 - adjusted_energy_loss,
            'mechanism': rejection_type,
            'section_damping': resonance.get('damping', 'unknown'),
            'explanation': f"Section R={r_value} dissipates {adjusted_energy_loss:.0%} of proposal energy"
        }
    
    def evaluate_proposal(self, proposal: Dict) -> Dict:
        """
        Evaluate a proposal against section-specific rules.
        
        Args:
            proposal: Proposal dictionary with section, type, reviews, etc.
        
        Returns:
            Evaluation result with approval status and metrics
        """
        section = proposal.get('section', '')
        proposal_type = proposal.get('type', 'addition')
        reviews = proposal.get('reviews', [])
        
        # Get threshold and rules
        threshold, rules = self.calculate_approval_threshold(section, proposal_type)
        
        # Calculate approval score
        approval_weights = rules.get('approval_weight', {
            'arbiter': 0.4,
            'reviewer': 0.3,
            'community': 0.3
        })
        
        total_score = 0
        total_weight = 0
        
        for review in reviews:
            role = review.get('role', 'reviewer')
            approval = review.get('approval', 0)  # 0-1 scale
            weight = approval_weights.get(role, approval_weights['reviewer'])
            
            total_score += approval * weight
            total_weight += weight
        
        # Normalize score
        approval_score = total_score / total_weight if total_weight > 0 else 0
        
        # Check other requirements
        enough_reviewers = len(reviews) >= rules.get('min_reviewers', 2)
        review_period_met = self._check_review_period(proposal, rules)
        quorum_met = True  # Implement actual quorum check if needed
        
        if rules.get('requires_quorum'):
            # Check if enough participants are active
            quorum_met = self._check_quorum(proposal)
        
        # Final approval decision
        approved = (
            approval_score >= threshold and
            enough_reviewers and
            review_period_met and
            quorum_met
        )
        
        return {
            'approved': approved,
            'approval_score': approval_score,
            'threshold': threshold,
            'enough_reviewers': enough_reviewers,
            'review_period_met': review_period_met,
            'quorum_met': quorum_met,
            'rules_applied': rules,
            'resonance': self.get_resonance_characteristics(section)
        }
    
    def _check_review_period(self, proposal: Dict, rules: Dict) -> bool:
        """Check if review period has elapsed."""
        if rules.get('fast_track'):
            return True
            
        submitted = proposal.get('submitted_at')
        if not submitted:
            return False
            
        # Parse timestamp
        if isinstance(submitted, str):
            submitted_dt = datetime.fromisoformat(submitted.replace('Z', '+00:00'))
        else:
            submitted_dt = submitted
            
        review_days = rules.get('review_period_days', 3)
        elapsed = datetime.now() - submitted_dt
        
        return elapsed >= timedelta(days=review_days)
    
    def _check_quorum(self, proposal: Dict) -> bool:
        """Check if enough participants are available for quorum."""
        # Simplified - would check LCT availability scores in real implementation
        return len(proposal.get('reviews', [])) >= 3
    
    def generate_proposal_template(self, section_path: str, 
                                  proposal_type: str = "addition") -> Dict:
        """
        Generate a proposal template with section-specific requirements.
        """
        threshold, rules = self.calculate_approval_threshold(section_path, proposal_type)
        resonance = self.get_resonance_characteristics(section_path)
        
        return {
            "proposal_id": None,  # To be assigned
            "section": section_path,
            "type": proposal_type,
            "title": "",
            "description": "",
            "content": "",
            "submitted_at": datetime.now().isoformat(),
            "requirements": {
                "approval_threshold": threshold,
                "min_reviewers": rules['min_reviewers'],
                "review_period_days": rules['review_period_days'],
                "token_cost": rules['token_cost'],
                "special_rules": rules.get('special_rules', [])
            },
            "resonance_characteristics": resonance,
            "status": "pending",
            "reviews": []
        }
    
    def validate_proposal_submission(self, proposal: Dict, 
                                    submitter_tokens: int) -> Tuple[bool, str]:
        """
        Validate if a proposal can be submitted.
        
        Args:
            proposal: Proposal to validate
            submitter_tokens: Available tokens of submitter
        
        Returns:
            Tuple of (is_valid, reason)
        """
        section = proposal.get('section', '')
        proposal_type = proposal.get('type', 'addition')
        
        # Get rules
        _, rules = self.calculate_approval_threshold(section, proposal_type)
        token_cost = rules.get('token_cost', 50)
        
        # Check token balance
        if submitter_tokens < token_cost:
            return False, f"Insufficient tokens: need {token_cost}, have {submitter_tokens}"
        
        # Check required fields
        required = ['title', 'description', 'content', 'section']
        for field in required:
            if not proposal.get(field):
                return False, f"Missing required field: {field}"
        
        # Check section exists
        section_path = self.project_root / "whitepaper" / "sections" / section
        if not section_path.exists():
            return False, f"Section does not exist: {section}"
        
        # Check special rules
        special_rules = rules.get('special_rules', [])
        if 'requires_philosophical_review' in special_rules:
            if 'philosophical_justification' not in proposal:
                return False, "This section requires philosophical justification"
        
        if 'requires_mathematical_validation' in special_rules:
            if 'mathematical_proof' not in proposal:
                return False, "This section requires mathematical validation"
        
        return True, "Valid"
    
    def report_governance_status(self) -> str:
        """Generate a governance status report with full LRC analysis."""
        report_lines = [
            "# Section Governance Configuration",
            f"Generated: {datetime.now().isoformat()}",
            "",
            "## LRC Resonance Spectrum",
            "",
            "| Section | L | C | R | Q | Damping | Frequency | Rejection |",
            "|---------|---|---|---|---|---------|-----------|-----------|"
        ]
        
        # Group sections by resonance
        resonance_groups = {}
        for section_name, section_rules in self.rules['sections'].items():
            resonance = section_rules.get('resonance', 'medium')
            if resonance not in resonance_groups:
                resonance_groups[resonance] = []
            resonance_groups[resonance].append({
                'name': section_name,
                'threshold': section_rules.get('change_threshold'),
                'description': section_rules.get('description', '')
            })
        
        # Sort resonance levels
        resonance_order = ['maximum', 'very-high', 'high', 'medium-high', 
                          'medium', 'low', 'very-low']
        
        for resonance in resonance_order:
            if resonance in resonance_groups:
                band = self.rules['resonance_bands'].get(resonance, {})
                
                for section in resonance_groups[resonance]:
                    report_lines.append(
                        f"| {section['name'][:20]} | "
                        f"{band.get('L', 'N/A')} | "
                        f"{band.get('C', 'N/A')} | "
                        f"{band.get('R', 'N/A')} | "
                        f"{band.get('Q_factor', 'N/A'):.2f} | "
                        f"{band.get('damping', 'N/A')[:8]} | "
                        f"{band.get('frequency', 'N/A')} | "
                        f"{band.get('rejection_rate', 'N/A'):.0%} |"
                    )
        
        report_lines.extend([
            "",
            "## Energy Dissipation Mechanisms",
            ""
        ])
        
        for mechanism_name, mechanism in self.rules.get('rejection_mechanisms', {}).items():
            report_lines.extend([
                f"### {mechanism_name.replace('_', ' ').title()}",
                f"- {mechanism['description']}",
                f"- Energy Loss: {mechanism['energy_loss']:.0%}",
                f"- Token Return: {mechanism['token_return']:.0%}",
                ""
            ])
        
        report_lines.extend([
            "## System Stability",
            "",
            "The R component ensures system stability by dissipating energy from bad proposals.",
            "Without R, the system would oscillate indefinitely or accumulate noise.",
            "",
            "- **Critical Damping (ζ=1)**: R = 2√(L/C) - smooth approach to equilibrium",
            "- **Overdamped (ζ>1)**: High R - slow response, high rejection",
            "- **Underdamped (ζ<1)**: Low R - oscillatory response, allows experimentation",
            ""
        ])
        
        return "\n".join(report_lines)


if __name__ == "__main__":
    # Test the governance system
    import sys
    
    project_root = Path(__file__).parent.parent.parent
    gov = SectionGovernance(project_root)
    
    if len(sys.argv) > 1:
        section = sys.argv[1]
        rules = gov.get_section_rules(section)
        print(f"Rules for {section}:")
        print(json.dumps(rules, indent=2))
        
        resonance = gov.get_resonance_characteristics(section)
        print(f"\nResonance characteristics:")
        print(json.dumps(resonance, indent=2))
    else:
        # Generate status report
        print(gov.report_governance_status())